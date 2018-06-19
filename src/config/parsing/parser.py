#!/usr/bin/env python
"""
    This script parses path information of Bioinformatics tools to a json Input File

                              Script Options

    -i      "The input Tool info text file "                                  (Required)
    -o      "json Input file to which path information is written into"       (Required)
"""

import json
import argparse
import sys
import string
import src.config.util.util as utility
from src.config.util.log import ProjectLogger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action='append', required=True,
                        help="The input configuration files. Multiple entries of this flag are allowed"
                        )
    parser.add_argument("--jsonTemplate", required=True,
                        help='The json template file to be filled in with data from the input files'
                        )
    parser.add_argument('--jobID', type=str, help='The job ID', default='NA', required=False)
    return parser.parse_args()


class Parser:
    def __init__(self, job_id="NA"):
        # Initialize the project logger
        self.project_logger = ProjectLogger(job_id, "parsing.parser.Parser")
        self.job_id = job_id

    @staticmethod
    def read_input_file(file_path):
        with open(file_path, "r") as F:
            return F.read().splitlines()

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

    def clean_input_file(self, input_lines):
        """
         Takes in a list of input lines, and removes any blank and comment lines (Lines beginning with '#')
        """
        # Remove all blank lines
        non_empty_lines = list(filter(None, input_lines))
        return self.remove_comments(non_empty_lines)

    def create_key_value_pairs(self, input_lines):
        """
        Turns a list of lines with keys and values separated by a '=' into pairs of (key, value) tuples
        """
        key_value_pairs = []

        for line in input_lines:
            if "=" not in line:
                sys.exit("ERROR: no equals sign present in the following line: " + line)
            else:
                # Get the position of the first equals sign (assumes there is no '=' in the key name)
                split_pos = line.index("=")
                key = line[0:split_pos]
                value = line[split_pos + 1:]
                key_value_pairs.append((key, value))
        return key_value_pairs

    def validate_key_value_pairs(self, key_value_pairs):
        """
         Takes in a list of (Key, Value) tuples, and confirms that they are valid (or throws an error)

         Checks performed:
            1. Verifies that all Keys have an associated Value
            2. Verifies that the Value is enclosed in double quotes
            3. Verifies that no special characters are present in the Values
            4. Verifies that no Key is present more than once
        """
        # List comprehensions to create individual list for Keys and Values
        keys = list(map(lambda x: x[0], key_value_pairs))
        values = list(map(lambda x: x[1], key_value_pairs))

        for key, value in list(zip(keys, values)):
            # Verifies if all the Tools have a corresponding Path to it
            if value == '':
                sys.exit("ERROR: Tool " + key + " is missing the path to the executable")

            # Check to see if the Paths are enclosed in double quotes in the Tool info text file
            elif value[0] != '"' or value[-1] != '"':
                sys.exit("ERROR: Executable for " + key + " is not enclosed in double quotes")

            # Check to verify if special characters are present in the Paths information
            for specialChar in string.punctuation.replace('/', '').replace('"', '').replace('-', '').replace('.', ''):
                if specialChar in value:
                    sys.exit("ERROR: Executable for " + key + ' has an invalid special character: "' + specialChar + '"')

            # Check whether any key is present multiple times
            if keys.count(key) > 1:
                sys.exit("ERROR: Tool " + key + " is listed twice in the tools list text file")

    def insert_values_into_dict(self, starting_dict, key_value_tuple):
        """
         Takes an initial dictionary and a list of (Key, Value) tuples.
           The Values in the tuple list are placed in the initial dictionary by Key.
        """
        output_dict = starting_dict.copy()

        # For each key-value pair that will be substituted into the dictionary
        for original_key, original_value in key_value_tuple:
            # Loop through each key, and pattern match to see if the original_key is found in this starting_dict.key
            for dict_key in starting_dict.keys():
                #  original_key matches the last section of this starting_dict.key
                #    keys are structured Major.Minor.KeyName, and we are trying to match against the KeyName
                if original_key == dict_key.split('.')[-1]:
                    # Trim quote marks off of original value
                    trimmed_value = original_value.replace('"', '')
                    output_dict[dict_key] = trimmed_value
        return output_dict

    def fill_in_json_template(self, input_file_list, json_template_file):
        """
         Takes in a list of input files and the location of the json file template, and writes an output file
           that contains the template's keys filled in with values from the input files

           The original template file will will be replaced with the filled-in version
        """
        # Read in the information from the json template file as a Python Dictionary
        #   The values of the template dictionary are filled in as input files are processed
        template_dict = utility.read_json_file(json_template_file, self.project_logger,
                                               json_not_found_error_code="placeholder_error_code",
                                               json_bad_format_error_code="placeholder_error_code"
                                               )
        # This loop iterates through all the input files and parses the path information
        for input_file in input_file_list:
            # Read in and clean the input file
            raw_input_lines = self.read_input_file(input_file)
            input_lines = self.clean_input_file(raw_input_lines)

            # Turn input lines into Tuples of Key-Value pairs
            key_value_tuples = self.create_key_value_pairs(input_lines)
            # Validate the key-value entries (Returns nothing; only possible outputs are error messages)
            self.validate_key_value_pairs(key_value_tuples)

            # Update the values in the template dictionary
            template_dict = self.insert_values_into_dict(template_dict, key_value_tuples)

        # Write the python dictionary out as a JSON file in the same location as the original template
        with open(json_template_file, "r+") as updated_json:
            json.dump(template_dict, updated_json, indent=4)


def main():
    args = parse_args()

    # List to hold all the input files
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

