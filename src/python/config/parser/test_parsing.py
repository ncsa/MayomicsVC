#!/usr/bin/env python3

import unittest
from config.parser.parsing import Parser


class TestParser(unittest.TestCase):

    # Create an instance of the Parser class
    parser_inst = Parser(job_id="NA")

    # Turn the project logger off during UnitTesting, so the end user is not confused by error messages
    #   (Some tests are designed to fail, so they will log "ERROR" messages that are expected)
    parser_inst.project_logger.logger.disabled = True

    def test_remove_comments(self):
        # Should remove comment lines
        input_lines = ["# Comment line", "      # Whitespace with comment", 'Key="Value"']
        filtered_lines = Parser.remove_comments(input_lines)
        self.assertEqual(filtered_lines, ['Key="Value"'])

    def test_clean_input_file(self):
        # Should remove blank and comment lines
        input_lines = ["", "", "# Comment line", 'Key="Value"']
        filtered_lines = Parser.clean_input_file(input_lines)
        self.assertEqual(filtered_lines, ['Key="Value"'])

    def test_create_key_value_pairs(self):
        # Note: the second test case purposefully has an '=' in the value (the parser only assumes the key has no '=')
        input_lines = ['Key1="Value1"', 'Key2="Value=2"']

        expected_output = [('Key1', '"Value1"'), ('Key2', '"Value=2"')]
        self.assertEqual(expected_output,
                         self.parser_inst.create_key_value_pairs(input_lines, "test_create_key_value_pairs")
                         )

    def test_validate_key_value_pairs_pass(self):
        '''
        This test has no assert. The method being tested returns nothing, but throws errors if anything fails

        This test should pass if the validate function can be called without throwing an error
        '''
        valid_tuple = [("keyA", '"valueA"')]
        self.parser_inst.validate_key_value_pairs(valid_tuple, file_path="dummy_file_path")

    def test_validate_key_value_pairs_fail_empty_value(self):
        no_value_tuple = [("keyA", "")]
        with self.assertRaises(SystemExit):
            self.parser_inst.validate_key_value_pairs(no_value_tuple, file_path="dummy_file_path")

    def test_validate_key_value_pairs_pass_empty_optional_key(self):
        # InputRead2 is a key that is allowed to be empty (see src/config/util/special_keys.py)
        nullable_key_empty_value = [("DebugMode", "")]
        self.parser_inst.validate_key_value_pairs(nullable_key_empty_value, file_path="dummy_file_path")

    def test_validate_key_value_pairs_fail_empty_non_optional_key(self):
        # InputRead1 is a key that is not allowed to be empty (it must have a value)
        key_empty_value = [("InputRead1", "")]
        with self.assertRaises(SystemExit):
            self.parser_inst.validate_key_value_pairs(key_empty_value, file_path="dummy_file_path")

    def test_validate_key_value_pairs_fail_no_quotes(self):
        no_value_tuple = [("keyA", 'Value without quotes')]
        with self.assertRaises(SystemExit):
            self.parser_inst.validate_key_value_pairs(no_value_tuple, file_path="dummy_file_path")

    def test_validate_key_value_pairs_fail_special_characters(self):
        no_value_tuple = [("keyA", '!@#$%&&^%(*&^%s')]
        with self.assertRaises(SystemExit):
            self.parser_inst.validate_key_value_pairs(no_value_tuple, file_path="dummy_file_path")

    def test_validate_key_value_pairs_fail_duplicate_keys(self):
        no_value_tuple = [("duplicateKey", 'valueA'), ("duplicateKey", "valueB")]
        with self.assertRaises(SystemExit):
            self.parser_inst.validate_key_value_pairs(no_value_tuple, file_path="dummy_file_path")

    def test_insert_values_into_dict(self):
        original_dict = {'major.minor.A': "init_A_value",
                         'major.minor.B': "init_B_value",
                         'major.minor.C': "init_C_value"
                         }
        key_value_tuples = [('A', '"final_A_value"'), ("B", '"final_B_value"')]

        substituted_dict = self.parser_inst.insert_values_into_dict(original_dict,
                                                                    key_value_tuples,
                                                               "test_insert_values_into_dict"
                                                                    )

        # The final dictionary should have new values for A and B, which C's value unchanged
        expected_dict = {'major.minor.A': "final_A_value",
                         'major.minor.B': "final_B_value",
                         'major.minor.C': "init_C_value"
                         }
        self.assertEqual(expected_dict, substituted_dict)

    def test_combine_input_read_arrays_paired_end_both(self):
        key_value_tuples = [("PairedEnd", '"true"'),
                            ("NormalInputRead1", '"readL1.fq,readL2.fq,readL3.fq"'),
                            ("NormalInputRead2", '"readR1.fq,readR2.fq,readR3.fq"')
                            ]
        expected_paired_end_value = [["readL1.fq", "readR1.fq"], ["readL2.fq", "readR2.fq"], ["readL3.fq", "readR3.fq"]]

        actual_paired_end_value = self.parser_inst.combine_input_read_arrays(key_value_tuples,
                                                                             "NormalInputRead1",
                                                                             "NormalInputRead2"
                                                                             )
        self.assertEqual(expected_paired_end_value, actual_paired_end_value)

    def test_combine_input_read_arrays_paired_end_one(self):
        key_value_tuples = [("PairedEnd", '"true"'),
                            ("NormalInputRead1", '"readL1.fq,readL2.fq,readL3.fq"'),
                            ("NormalInputRead2", '""')
                            ]
        with self.assertRaises(SystemExit):
            # Should fail, as paired end is true but only one read set is provided
            self.parser_inst.combine_input_read_arrays(key_value_tuples, "NormalInputRead1", "NormalInputRead2")

    def test_combine_input_read_arrays_paired_end_unequal_lists(self):
        key_value_tuples = [("PairedEnd", '"true"'),
                            ("NormalInputRead1", '"readL1.fq,readL2.fq,readL3.fq"'),
                            ("NormalInputRead2", '"readR1.fq"')
                            ]
        with self.assertRaises(SystemExit):
            # Should fail, as paired end is true but only one read set is provided
            self.parser_inst.combine_input_read_arrays(key_value_tuples, "NormalInputRead1", "NormalInputRead2")

    def test_combine_input_read_arrays_single_end_both(self):
        key_value_tuples = [("PairedEnd", '"false"'),
                            ("NormalInputRead1", '"readL1.fq,readL2.fq,readL3.fq"'),
                            ("NormalInputRead2", '"readR1.fq,readR2.fq,readR3.fq"')
                            ]
        expected_paired_end_value = [["readL1.fq"], ["readL2.fq"], ["readL3.fq"],
                                     ["readR1.fq"], ["readR2.fq"], ["readR3.fq"]]

        actual_paired_end_value = self.parser_inst.combine_input_read_arrays(key_value_tuples,
                                                                             "NormalInputRead1",
                                                                             "NormalInputRead2"
                                                                             )
        self.assertEqual(expected_paired_end_value, actual_paired_end_value)

    def test_combine_input_read_arrays_single_end_one(self):
        key_value_tuples = [("PairedEnd", '"false"'),
                            ("NormalInputRead1", '"readL1.fq,readL2.fq,readL3.fq"'),
                            ("NormalInputRead2", '""')
                            ]
        expected_paired_end_value = [["readL1.fq"], ["readL2.fq"], ["readL3.fq"]]

        actual_paired_end_value = self.parser_inst.combine_input_read_arrays(key_value_tuples,
                                                                             "NormalInputRead1",
                                                                             "NormalInputRead2"
                                                                             )
        self.assertEqual(expected_paired_end_value, actual_paired_end_value)



