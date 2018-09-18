#!/usr/bin/env python
"""
    This script parses the values from a flat 'key="value"' formatted file into a JSON template file by key

                              Script Options

    -i                  "The input Tool info text file"                                  (Required)
    --jsonTemplate      "json Input file to which path information is written into"      (Required)

    --jobID             "The job ID (defaults to 'NA')"                                  (Optional)
    -d                  "Enables debug mode"                                             (Optional)
"""

import json
import argparse
import logging
import sys

from ..util.util import read_json_file
from ..util.log import ProjectLogger
from ..util.special_keys import OPTIONAL_KEYS

"""
Exit code Rules:

1. Exit codes in this module are only given when an error has occurred, so they will all start with 'E.'
2. The letters 'par.' because they are coming from the parser component of the code
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
"""


def parse_args():
    """
    By default, argparse treats all arguments that begin with '-' or '--' as optional in the help menu
      (preferring to have required arguments be positional).

    To get around this, we must define a required group to contain the required arguments
      This will cause the help menu to be displayed correctly
    """
    parser = argparse.ArgumentParser()

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("-i", action='append', required=True, metavar='',
                                help="The input configuration files (Multiple entries of this flag are allowed)"
                                )
    required_group.add_argument("--jsonTemplate", required=True, metavar='',
                                help='The json template file that is filled in with data from the input files'
                                )
    required_group.add_argument("-o", required=True, metavar='',
                                help='The location of the output file'
                                )
    # Truly optional argument
    parser.add_argument('--jobID', type=str, metavar='', help='The job ID', default='NA', required=False)
    # Debug mode is on when the flag is present and is false by default
    parser.add_argument("-d", action="store_true", help="Turns on debug mode", default=False, required=False)
    return parser.parse_args()


class Parser:
    def __init__(self, job_id="NA", debug_mode=False):
        # Initialize the project logger
        if debug_mode:
            self.project_logger = ProjectLogger(job_id, "parsing.parser.Parser", logging.DEBUG)
        else:
            self.project_logger = ProjectLogger(job_id, "parsing.parser.Parser")
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

    def handle_special_keys(self, config_key, full_dict_key, output_dictionary, trimmed_value, file_path):
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
                    "REF key from input file '" + file_path +
                    "' appears to not have a valid value: '" + trimmed_value + "has no '.fa' or '.fasta' extension"
                )
                sys.exit(1)
            else:
                base_name = trimmed_value[:extension_start_index]
                output_dictionary[full_dict_key + "Dict"] = base_name + ".dict"

        elif config_key == "DBSNP":
            output_dictionary[str(full_dict_key) + "Idx"] = str(trimmed_value) + ".idx"

    def insert_values_into_dict(self, starting_dict, key_value_tuple, file_path):
        """
         Takes an initial dictionary and a list of (Key, Value) tuples.
           The Values in the tuple list are placed in the initial dictionary by Key.
        """
        output_dict = starting_dict.copy()

        # For each key-value pair that will be substituted into the dictionary
        for config_key, config_value in key_value_tuple:
            # Switch signaling whether the key from the config file was found in the json template
            config_key_was_present = False

            # Loop through each key, and pattern match to see if the config_key is found in this starting_dict.key
            for dict_key in starting_dict.keys():
                #  config_key matches the last section of this starting_dict.key
                #    keys are structured Major.Minor.KeyName, and we are matching against the KeyName
                if config_key == dict_key.split('.')[-1]:
                    config_key_was_present = True

                    # Remove quote marks from the ends of the original value; assumes that the value was wrapped in
                    #   quote marks (if it got past the validate_key_value_pairs method, this is guaranteed)
                    trimmed_value = config_value[1:-1]

                    output_dict[dict_key] = trimmed_value

                    # Handle special keys that need additional json keys added for each config key (such as REF, DBSNP)
                    self.handle_special_keys(config_key, dict_key, output_dict, trimmed_value, file_path)

            # Log a warning message if a key was in the config file but not in the template
            if not config_key_was_present:
                self.project_logger.log_warning(
                    "Key '" + config_key + "' in config file '" + file_path +
                    "' had no corresponding key in the JSON template; this key-value pair was ignored"
                )
        return output_dict

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
        # This loop iterates through all the input files and parses the path information
        for input_file in input_file_list:
            # Read in and clean the input file
            raw_input_lines = self.read_input_file(input_file)
            input_lines = self.clean_input_file(raw_input_lines)

            # Turn input lines into Tuples of Key-Value pairs
            key_value_tuples = self.create_key_value_pairs(input_lines, file_path=input_file)
            # Validate the key-value entries (Returns nothing; only possible outputs are error messages)
            self.validate_key_value_pairs(key_value_tuples, file_path=input_file)

            # Update the values in the template dictionary
            template_dict = self.insert_values_into_dict(template_dict, key_value_tuples, file_path=input_file)

        # Write the python dictionary out as a JSON file in the output file location
        with open(output_file, "w") as updated_json:
            json.dump(template_dict, updated_json, indent=4, sort_keys=True)

        # Write a success message to the log
        self.project_logger.log_info(
                'Configuration file parsing finished successfully with ' + str(self.project_logger.warnings_issued) +
                ' warning(s) issued'
        )


def main():
    args = parse_args()

    # List of all the input files
    input_file_list = args.i

    # Instantiation of the Parser class
    if args.jobID is None:
        k_v_parser = Parser(debug_mode=args.d)
    else:
        k_v_parser = Parser(args.jobID, debug_mode=args.d)

    # Fill in the json template file with values from the Key="Value" formatted input files
    k_v_parser.fill_in_json_template(input_file_list, args.jsonTemplate, args.o)


if __name__ == '__main__':
    main()

