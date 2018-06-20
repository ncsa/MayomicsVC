#!/usr/bin/env python
"""
    This script parses the values from a flat 'key="value"' formatted file into a JSON template file by key

                              Script Options

    -i                  "The input Tool info text file"                                  (Required)
    --jsonTemplate      "json Input file to which path information is written into"      (Required)

    --jobID             "The job ID (defaults to 'NA')"                                  (Optional)
"""

import json
import argparse
import sys
import string
import src.config.util.util as utility
from src.config.util.log import ProjectLogger

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
E.par.SpC.1 = A non-comment line in a config file had an invalid special character
E.par.Key.1 = A key is present multiple times in a config file

E.par.JSN.1 = An input JSON file could not be found
E.par.JSN.2 = An input JSON file was not formatted properly
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
                                help='The json template file to be filled in with data from the input files'
                                )
    # Truly optional argument
    parser.add_argument('--jobID', type=str, metavar='', help='The job ID', default='NA', required=False)
    return parser.parse_args()


class Parser:
    def __init__(self, job_id="NA"):
        # Initialize the project logger
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

         Checks performed:
            1. Verifies that all Keys have an associated Value
            2. Verifies that the Value is enclosed in double quotes
            3. Verifies that no special characters are present in the Values
            4. Verifies that no Key is present more than once
        """
        # List comprehensions to create complete list of keys
        keys_list = [k for k, v in key_value_pairs]

        for key, value in key_value_pairs:
            # Verifies if all the Tools have a corresponding Path to it
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
            # Check to verify if special characters are present in the Paths information
            for specialChar in string.punctuation.replace('/', '').replace('"', '').replace('-', '').replace('.', ''):
                if specialChar in value:
                    self.project_logger.log_error(
                        "E.par.SpC.1",
                        "Invalid special character '" + specialChar + "' found in value '" + value +
                        "' of key '" + key + "' in input file '" + file_path + "'"
                    )
                    sys.exit(1)
            # Check whether any key is present multiple times
            if keys_list.count(key) > 1:
                self.project_logger.log_error(
                    "E.par.Key.1", "Key '" + key + "' is present more than once in input file '" + file_path + "'"
                )
                sys.exit(1)

    def insert_values_into_dict(self, starting_dict, key_value_tuple, file_path):
        """
         Takes an initial dictionary and a list of (Key, Value) tuples.
           The Values in the tuple list are placed in the initial dictionary by Key.
        """
        output_dict = starting_dict.copy()

        # For each key-value pair that will be substituted into the dictionary
        for original_key, original_value in key_value_tuple:
            # Switch signaling whether the key from the config file was found in the json template
            config_key_was_present = False

            # Loop through each key, and pattern match to see if the original_key is found in this starting_dict.key
            for dict_key in starting_dict.keys():
                #  original_key matches the last section of this starting_dict.key
                #    keys are structured Major.Minor.KeyName, and we are matching against the KeyName
                if original_key == dict_key.split('.')[-1]:
                    config_key_was_present = True
                    # Trim quote marks off of original value
                    trimmed_value = original_value.replace('"', '')
                    output_dict[dict_key] = trimmed_value

            # Log a warning message if a key was in the config file but not in the template
            if not config_key_was_present:
                self.project_logger.log_warning(
                    "Key '" + original_key + "' in config file '" + file_path +
                    "' had no corresponding key in the JSON template; this key-value pair was ignored"
                )

        return output_dict

    def fill_in_json_template(self, input_file_list, json_template_file):
        """
         Takes in a list of input files and the location of the json file template, and writes an output file
           that contains the template's keys filled in with values from the input files

         The original template file will be replaced with the filled-in version
        """
        # Read in the information from the json template file as a Python Dictionary
        #   The values of the template dictionary are filled in as input files are processed
        template_dict = utility.read_json_file(json_template_file, self.project_logger,
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

        # Write the python dictionary out as a JSON file in the same location as the original template
        with open(json_template_file, "r+") as updated_json:
            json.dump(template_dict, updated_json, indent=4)

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
        k_v_parser = Parser()
    else:
        k_v_parser = Parser(args.jobID)

    # Fill in the json template file with values from the Key="Value" formatted input files
    k_v_parser.fill_in_json_template(input_file_list, args.jsonTemplate)


if __name__ == '__main__':
    main()

