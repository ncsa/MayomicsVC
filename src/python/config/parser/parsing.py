#!/usr/bin/env python3

import sys
from util.util import read_json_file
from util.log import ProjectLogger
from config.util.special_keys import OPTIONAL_KEYS
import json
import logging
from typing import Dict

"""
Exit code Rules:

1. Exit codes in this module are only given when an error has occurred, so they will all start with 'E.'
2. The letters 'par.' because they are coming from the parser_inst component of the code
3. A three letter code that hints at what the problem was
4. A number that can help to differentiate similar error codes


Error Code List
========================================================================================================================
E.par.Fil.1 = An input file could not be found
E.par.NEq.1 = A non-comment line in a config file had no equals sign 
E.par.NVa.1 = A non-comment line in a config file had no value specified
E.par.NQt.1 = A non-comment line in a config file had a value that was not enclosed in quotes
E.par.WhS.1 = A non-comment line in a config file had a value with only whitespace in between the quote marks
E.par.SpC.1 = A non-comment line in a config file had an invalid special character
E.par.Key.1 = A key is present multiple times in a config file

E.par.JSN.1 = An input JSON file could not be found
E.par.JSN.2 = An input JSON file was not formatted properly

E.par.REF.1 = The config file REF key points to a value that does not have an '.fa' or '.fasta' extension

E.par.InR.1 = Either the PairedEnd, InputRead1, or InputRead2 key was not present in the config file
E.par.InR.2 = The InputRead1 and InputRead2 lists had different lengths
E.par.InR.3 = PairedEnd was true but the InputRead2 lists was empty
"""


class Parser:
    def __init__(self, job_id="NA", debug_mode=False):
        # Initialize the project logger
        if debug_mode:
            self.project_logger = ProjectLogger(job_id, "parser.parsing.Parser", logging.DEBUG)
        else:
            self.project_logger = ProjectLogger(job_id, "parser.parsing.Parser")
        self.job_id = job_id

    def read_input_file(self, file_path):
        try:
            with open(file_path, "r") as F:
                return F.read().splitlines()
        except FileNotFoundError:
            self.project_logger.log_error("E.par.Fil.1", 'Input file "' + str(file_path) + '" could not be found')
            sys.exit(1)

    @staticmethod
    def remove_comments(input_lines):
        """
         Remove any comment lines from the list of lines

           A comment line is a line where the first non-whitespace character is a '#'
        """
        filtered_lines = []
        for line in input_lines:
            if line != "":
                # If the first non-space character is a '#', exclude it
                if line.strip()[0] != '#':
                    filtered_lines.append(line)
        return filtered_lines

    @staticmethod
    def clean_input_file(input_lines):
        """
         Takes in a list of input lines, and removes any blank and comment lines (lines beginning with '#')
        """
        # Remove all blank lines
        non_empty_lines = list(filter(None, input_lines))
        return Parser.remove_comments(non_empty_lines)

    def create_key_value_pairs(self, input_lines, file_path):
        """
        Turns a list of lines with keys and values separated by a '=' into pairs of (key, value) tuples
        """
        key_value_pairs = []

        for line in input_lines:
            if "=" not in line:
                self.project_logger.log_error(
                    "E.par.NEq.1",
                    "No equals sign present in line '" + line + "' from input file '" + file_path + "'"
                )
                sys.exit(1)
            else:
                # Get the position of the first equals sign (assumes there is no '=' in the key name)
                split_pos = line.index("=")
                key = line[0:split_pos]
                value = line[split_pos + 1:]
                key_value_pairs.append((key, value))
        return key_value_pairs

    def validate_key_value_pairs(self, key_value_pairs, file_path):
        """
         Takes in a list of (Key, Value) tuples, and confirms that they are valid (or throws an error)

        Keys that are allowed to be empty are not checked (see src/config/util/special_keys.py)

         Checks performed:
            1. Verifies that all Keys have an associated Value
            2. Verifies that the Value is enclosed in double quotes
            3. Verifies that no special characters are present in the Values
            4. Verifies that no Key is present more than once
            5. Verifies that no value is whitespace only
        """
        # List comprehensions to create complete list of keys
        keys_list = [k for k, v in key_value_pairs]

        for key, value in key_value_pairs:
            if key.lower() in OPTIONAL_KEYS and (value == "" or value == '""'):
                # These keys are allowed to have empty values, do not perform checking (simply write a debug message)
                self.project_logger.log_debug(
                    "The key '" + key + "' had an empty value; since its value is optional, no error was thrown"
                )
            else:
                # Check that the value is not empty
                if value == '':
                    self.project_logger.log_error(
                        "E.par.NVa.1",
                        "No value present for key '" + key + "' in input file '" + file_path + "'"
                    )
                    sys.exit(1)
                # Check that the value is enclosed in double quotes
                elif value[0] != '"' or value[-1] != '"':
                    self.project_logger.log_error(
                        "E.par.NQt.1",
                        "No quotes around the value for key '" + key + "' in input file '" + file_path + "'"
                    )
                    sys.exit(1)
                # Check to see that non-whitespace are present between the quote marks
                #   value[1:-1] trims off the first and last chars and strip removes all whitespace chars from the ends
                elif value[1:-1].strip() == '':
                    self.project_logger.log_error(
                        "E.par.WhS.1",
                        "Only whitespace found in value '" + value + "' of key '" + key + "' in input file '" +
                        file_path + "'"
                    )
                    sys.exit(1)
                # Check if any special characters are present
                special_chars = "!#$%&()*;<>?@[]^`{|}~"
                for special_char in special_chars:
                    if special_char in value:
                        self.project_logger.log_error(
                            "E.par.SpC.1",
                            "Invalid special character '" + special_char + "' found in value '" + value +
                            "' of key '" + key + "' in input file '" + file_path + "'"
                        )
                        sys.exit(1)
                # Check whether any key is present multiple times
                if keys_list.count(key) > 1:
                    self.project_logger.log_error(
                        "E.par.Key.1", "Key '" + key + "' is present more than once in input file '" + file_path + "'"
                    )
                    sys.exit(1)
                else:
                    self.project_logger.log_debug("The key-value pair '" + key + "=" + value + "' is a valid pair")

    def handle_special_keys(self, config_key, full_dict_key, output_dictionary, trimmed_value):
        """
         Wrapper function designed to handle all special keys, i.e. configuration keys that need more than one
           corresponding value added to the JSON config file

         Designed to be extendable as the workflow changes

         Example:
            # In the config file, there is only one key
            DBSNP="dbsnp.vcf"

            # But in the JSON file, we need extra information:
            {
                "major.minor.DBSNP": "dbsnp.vcf"
                "major.minor.DBSNP_IDX": "dbsnp.vcf.idx"
                ...
            }

        :side-effect adds any additional keys derived from the input config_key to the output dictionary
        """
        def add_index_key_value(key_suffix, value_suffix):
            json_key = full_dict_key + key_suffix
            json_value = trimmed_value + "." + value_suffix
            output_dictionary[json_key] = json_value

        if config_key == "Ref":
            add_index_key_value("Amb", "amb")
            add_index_key_value("Ann", "ann")
            add_index_key_value("Bwt", "bwt")
            add_index_key_value("Pac", "pac")
            add_index_key_value("Sa", "sa")

            # dict files replace the '.fasta'/'.fa' extension with '.dict'
            # 'find' returns the index of the first char in the string searched for (or -1 if the string was not found)
            extension_start_index = trimmed_value.find(".fa")
            if extension_start_index == -1:
                self.project_logger.log_error(
                    "E.par.REF.1",
                    "REF key appears to not have a valid value: '" + trimmed_value +
                    " has no '.fa' or '.fasta' extension"
                )
                sys.exit(1)
            else:
                base_name = trimmed_value[:extension_start_index]
                output_dictionary[full_dict_key + "Dict"] = base_name + ".dict"

        elif config_key == "DBSNP":
            output_dictionary[str(full_dict_key) + "Idx"] = str(trimmed_value) + ".idx"

    def combine_input_read_arrays(self, key_value_tuples, read1_name, read2_name):
        """
        Takes in the key value tuples, and the names of two read variables, and combines their values into a 2D array
        (This is used for two pairs of reads: NormalInputReads and TumorInputReads)

        InputReads1 & 2 are each provided in the config file as an array of files

            InputRead1="L1,L2,L3, ..."
            InputRead2="R1,R2,R3, ..."

        These two lists will be combined into a single 2D-array variable called InputReads.

        How they are combined depends on the value of the PairedEnd flag

        If PairedEnd and data in both input variables
            InputRead = [ [L1, R1], [L2, R2], [L3, R3] ]

        If not PairedEnd and data in both input variables
            InputRead = [ [L1], [L2], [L3], [R1], [R2], [R3] ]

        If not PairedEnd and data only in InputRead1
            InputRead = [ [L1], [L2], [L3] ]

        Throw an error if:
            1. PairedEnd, InputRead1, or InputRead2 are not defined in the input tuple list
            1. PairedEnd is true and data only InputRead1 is set
            2. PairedEnd is true and the number of entries between the input variables are different

        :return: The 2-D array holding the values that will be the value of the InputReads key in the JSON file
        """
        try:
            # next iterates through a collection until it finds an item that meets a criterion; here it is used to
            #   find the tuples with the keys we want
            paired_end_tuple = next(pair for pair in key_value_tuples if pair[0] == "PairedEnd")
            input_read_1_tuple = next(pair for pair in key_value_tuples if pair[0] == read1_name)
            input_read_2_tuple = next(pair for pair in key_value_tuples if pair[0] == read2_name)
        except StopIteration:
            self.project_logger.log_error(
                "E.par.InR.1",
                "Either the PairedEnd, " + read1_name + ", or " + read2_name + " key was not present in the config file"
            )
            sys.exit(1)

        paired_end = True if paired_end_tuple[1].strip('"') == "true" else False

        # Split by comma (unless the string is empty, then return an empty list)
        input_read_1_value = input_read_1_tuple[1].strip('"')
        input_read_2_value = input_read_2_tuple[1].strip('"')

        input_read_1_array = [] if input_read_1_value == "" else input_read_1_value.split(",")
        input_read_2_array = [] if input_read_2_value == "" else input_read_2_value.split(",")

        if paired_end:
            if len(input_read_2_array) == 0:
                self.project_logger.log_error("E.par.InR.3", "PairedEnd was true but the InputRead2 lists was empty")
                sys.exit(1)
            elif len(input_read_1_array) != len(input_read_2_array):
                self.project_logger.log_error("E.par.InR.2", "The InputRead1 & InputRead2 lists had different lengths")
                sys.exit(1)

            # By default, zip returns tuples. The list comprehension turns them back into lists
            input_reads_2d_array = [list(x) for x in list(zip(input_read_1_array, input_read_2_array))]
        else:
            # Concatenate the two lists, then turn each entry into its own list
            input_reads_2d_array = [[x] for x in input_read_1_array + input_read_2_array]

        return input_reads_2d_array

    def insert_values_into_dict(self,
                                starting_dict,
                                key_value_tuple,
                                normal_input_reads_2d_array,
                                tumor_input_reads_2d_array=None
                                ):
        """
         Takes an initial dictionary and a list of (Key, Value) tuples.
           The Values in the tuple list are placed in the initial dictionary by Key.

        :param starting_dict = The dictionary derived from the original JSON template
        :param key_value_tuple = The complete list of keys and values (as tuples) that were in the config files
        :param normal_input_reads_2d_array = The 2D array that will be the value of the JSON "NormalInputReads" key
        :param tumor_input_reads_2d_array = 2D array that will be the value of the "TumorInputReads" key (if it is used)

        :return The dictionary with all of the values inserted
        """
        output_dict = starting_dict.copy()

        # Switches so internal if statements are not run repeatedly
        normal_input_reads_found = False
        tumor_input_reads_found = False

        # For each key-value pair that will be substituted into the dictionary
        for config_key, config_value in key_value_tuple:
            # Switch signaling whether the key from the config file was found in the json template
            config_key_was_present = False

            # Loop through each key, and pattern match to see if the config_key is found in this starting_dict.key
            for dict_key in starting_dict.keys():
                #  config_key matches the last section of this starting_dict.key
                #    keys are structured Major.Minor.KeyName, and we are matching against the KeyName
                dict_key_suffix = dict_key.split(".")[-1]

                if config_key == dict_key_suffix:
                    config_key_was_present = True

                    # Remove quote marks from the ends of the original value; assumes that the value was wrapped in
                    #   quote marks (if it got past the validate_key_value_pairs method, this is guaranteed)
                    trimmed_value = config_value[1:-1]

                    # Special case where the value should be split into an array
                    if dict_key_suffix == "PlatformUnit":
                        output_dict[dict_key] = trimmed_value.split(",")
                    else:
                        output_dict[dict_key] = trimmed_value

                    # Handle special keys that need additional json keys added for each config key (such as REF, DBSNP)
                    self.handle_special_keys(config_key, dict_key, output_dict, trimmed_value)

                # Add the special key, NormalInputReads, to the dictionary
                elif dict_key_suffix == "NormalInputReads" and not normal_input_reads_found:
                    output_dict[dict_key] = normal_input_reads_2d_array
                    # Flip the switch so that this conditional is not evaluated again
                    normal_input_reads_found = True

                elif dict_key_suffix == "TumorInputReads" and \
                        tumor_input_reads_2d_array is not None and \
                        not tumor_input_reads_found:
                    self.project_logger.log_debug(
                        "The TumorInputReads JSON key was paired with a 2D array created from the " +
                        "TumorInputRead1 and 2 variables"
                    )
                    output_dict[dict_key] = tumor_input_reads_2d_array
                    # Flip the switch so that this conditional is not evaluated again
                    tumor_input_reads_found = True

            # Log a warning message if a key was in the config file but not in the template
            exception_list = ["NormalInputRead1", "NormalInputRead2", "TumorInputRead1", "TumorInputRead2"]
            if not config_key_was_present and config_key not in exception_list:
                self.project_logger.log_warning(
                    "Key '" + config_key + "' had no corresponding key in the JSON template;" +
                    " this key-value pair was ignored"
                )
        return output_dict

    def is_TumorInputRead1_present(self, key_value_tuples) -> bool:
        """
        Checks whether the TumorInputRead1 key is present and if it has a value assigned to it

        :return: True if the variable has a value and False if it has no value or is not present
        """
        try:
            tumor_input_read_1_tuple = next(pair for pair in key_value_tuples if pair[0] == "TumorInputRead1")
            if tumor_input_read_1_tuple[1] in ("", '""'):
                self.project_logger.log_debug("The TumorInputRead1 key was found, and its value was empty")
                return False
            else:
                self.project_logger.log_debug("The TumorInputRead1 key was found, and its value was not empty")
                return True
        except StopIteration:
            self.project_logger.log_debug("The TumorInputRead1 key was not found in the config file")
            # The variable was not present
            return False

    def find_variables_in_JSON_not_in_config(self, all_config_tuples, JSON_dict: Dict[str, str]):
        """
        Warns the user if any JSON key was not present in the configuration tuple list (that are not in the exception
          list).

        Since there are some keys that are added to the JSON that are not part of the config file directly (see the
          "handle_special_keys" method), they are not sought for in the config keys list, and are put in the exception
          list

        :param all_config_tuples: list of all config file tuples
        :param JSON_dict: Python dictionary from the template JSON file
        """
        exception_list = ("RefAmb", "RefAnn", "RefBwt", "RefDict", "RefPac", "RefSa",
                          "DBSNPIdx",
                          "NormalInputReads", "TumorInputReads"
                          )

        # The informative part of the JSON key is the subsection after the last '.' character
        trimmed_JSON_keys = list(set([i.split(".")[-1] for i in JSON_dict.keys()]))
        config_keys = [i[0] for i in all_config_tuples]

        for json_key in trimmed_JSON_keys:
            if json_key not in exception_list and json_key not in config_keys:
                self.project_logger.log_warning("The '" + json_key + "' key in the JSON template did not have a " +
                                                "corresponding key in any of the config files; " +
                                                "this key was not filled in"
                                                )


    def fill_in_json_template(self, input_file_list, json_template_file, output_file):
        """
         Takes in a list of input files and the location of the json file template, and writes an output file
           that contains the template's keys filled in with values from the input files

         The original template files contents will be copied, filled-in, and saved in the new output file
        """
        # Read in the information from the json template file as a Python Dictionary
        #   The values of the template dictionary are filled in as input files are processed
        template_dict = read_json_file(json_template_file, self.project_logger,
                                       json_not_found_error_code="E.par.JSN.1",
                                       json_bad_format_error_code="E.par.JSN.2"
                                       )

        # Combine all of the input file key-value tuples into a single list
        all_key_value_tuples = []

        # This loop iterates through all the input files and parses the path information
        for input_file in input_file_list:
            # Read in and clean the input file
            raw_input_lines = self.read_input_file(input_file)
            input_lines = self.clean_input_file(raw_input_lines)

            # Turn input lines into Tuples of Key-Value pairs
            key_value_tuples = self.create_key_value_pairs(input_lines, file_path=input_file)
            # Validate the key-value entries (Returns nothing; only possible outputs are error messages)
            self.validate_key_value_pairs(key_value_tuples, file_path=input_file)

            # Add all of the pairs from this file into the overall list of pairs
            [all_key_value_tuples.append(pair) for pair in key_value_tuples]

        # Create the NormalInputReads value from the info in the config files (Always required)
        normal_input_reads_2d_array = self.combine_input_read_arrays(all_key_value_tuples,
                                                                     "NormalInputRead1",
                                                                     "NormalInputRead2"
                                                                     )

        # Optionally, create the TumorInputReads value from the info in the config files (will be a 2D array if
        #   TumorInputRead1 was not empty, and None otherwise)
        if self.is_TumorInputRead1_present(all_key_value_tuples):
            tumor_input_reads_2d_array = self.combine_input_read_arrays(all_key_value_tuples,
                                                                        "TumorInputRead1",
                                                                        "TumorInputRead2"
                                                                        )
            template_dict = self.insert_values_into_dict(template_dict,
                                                         all_key_value_tuples,
                                                         normal_input_reads_2d_array,
                                                         tumor_input_reads_2d_array
                                                         )
        else:
            # Update the values in the template dictionary
            template_dict = self.insert_values_into_dict(template_dict,
                                                         all_key_value_tuples,
                                                         normal_input_reads_2d_array
                                                         )

        # Send a warning if any JSON keys had no corresponding key in any of the config files
        self.find_variables_in_JSON_not_in_config(all_key_value_tuples, template_dict)

        # Write the python dictionary out as a JSON file in the output file location
        with open(output_file, "w") as updated_json:
            json.dump(template_dict, updated_json, indent=4, sort_keys=True)

        # Write a success message to the log
        self.project_logger.log_info(
                'Configuration file parser finished successfully with ' + str(self.project_logger.warnings_issued) +
                ' warning(s) issued'
        )
