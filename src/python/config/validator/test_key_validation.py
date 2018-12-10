#!/usr/bin/env python

import os
import sys
import unittest
from config.validator.key_validation import Validator

# This is the full path to the validation package
package_full_path = os.path.abspath(os.path.dirname(__file__))


class TestValidator(unittest.TestCase):

    # Create an instance of the validator (without passing in a job_id string, it defaults to "NA")
    validator = Validator()
    # Turn the project logger off during UnitTesting, so the end user is not confused by error messages
    #   (Some tests are designed to fail, so they will log "ERROR" messages that are expected)
    validator.project_logger.logger.disabled = True

    def test_trim_config_file_keys(self):
        full_config = {"major.minor.KeyName1": 1, "major.minor.KeyName2": 2}
        expected_trimmed_config = {"KeyName1": 1, "KeyName2": 2}

        actual_trimmed_config = self.validator.trim_config_file_keys(full_config)
        self.assertEqual(expected_trimmed_config, actual_trimmed_config)

    def test_check_key_file(self):
        # The test file is the key_validation.py source file
        test_file = package_full_path + "/key_validation.py"

        result = self.validator.check_key(key_name="key_to_file", key_value=test_file, key_type="File")
        self.assertTrue(result)

    def test_check_key_execfile(self):
        # The test file is the python executable used to call this script
        test_exec_file = sys.executable

        result = self.validator.check_key(key_name="key_to_execfile", key_value=test_exec_file, key_type="ExecFile")
        self.assertTrue(result)

    '''
    Convenience function to enable easy testing of different key types 
    
      Takes a test value and the expected type of the value (as understood by key_validation.Validator.check_keys()
        method)
      
      An optional parameter 'expected_to_fail' describes the type of test being conducted. By default, it expects the
      key check to return successfully, but if 'expected_to_fail' is True, the test is actually supposed to fail.
    '''
    def __check_key_tester(self, value, key_type, expected_to_fail=False):
        result = self.validator.check_key("test_key_name", value, key_type)

        if expected_to_fail:
            self.assertFalse(
                result,
                msg=value + " was considered a valid " + key_type + " even though it should not be"
            )
        else:
            self.assertTrue(result, msg=value + " is not considered a valid " + key_type + " type")

    def test_check_key_int(self):
        self.__check_key_tester(value="8", key_type="Integer")

    def test_check_key_int_failure(self):
        self.__check_key_tester(value="NotAnInteger", key_type="Integer", expected_to_fail=True)

    def test_check_key_bool(self):
        self.__check_key_tester(value="true", key_type="Boolean")
        self.__check_key_tester(value="false", key_type="Boolean")

    def test_check_key_bool_failure(self):
        self.__check_key_tester(value="True", key_type="Boolean", expected_to_fail=True)

    def test_check_key_decimal(self):
        self.__check_key_tester(value="3.14159", key_type="Decimal")
        self.__check_key_tester(value="10", key_type="Decimal")

    def test_check_key_decimal_failure(self):
        self.__check_key_tester(value="Nineteen Eighty-four", key_type="Decimal", expected_to_fail=True)

    def test_check_key_unrecognized_type(self):
        # An unrecognized key type should return false from the checker method
        result = self.validator.check_key("key_of_unknown_type", "KeyValue", key_type="UnacceptableType")
        self.assertFalse(result)

    def test_check_key_optional_type_with_empty_value(self):
        # DebugMode is listed as an optional key in src/config/util/special_keys.py, so it can be empty
        result = self.validator.check_key("DebugMode", "", key_type="DebugMode")
        self.assertTrue(result)


if __name__ == "__main__":
    unittest.main()
