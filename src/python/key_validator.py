#!/usr/bin/env python3

import argparse
import sys

from util.util import read_json_file
from config.validator.key_validation import Validator


def parse_args(args):
    """
    By default, argparse treats all arguments that begin with '-' or '--' as optional in the help menu
        (preferring to have required arguments be positional).

    To get around this, we must define a required group to contain the required arguments
        This will cause the help menu to be displayed correctly
    """
    parser = argparse.ArgumentParser(description='Validate values in Cromwell/WDL input file')

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('-i', metavar="", type=str, help='JSON configuration file to validate', required=True)
    required_group.add_argument('--KeyTypeFile',
                                type=str,
                                metavar="",
                                help='JSON file with typing info for keys',
                                required=True
                                )

    # Truly optional argument
    parser.add_argument('--jobID', type=str, metavar='', help='The job ID', default='NA', required=False)
    # Debug mode is on when the flag is present and is false by default
    parser.add_argument("-d", action="store_true", help="Turns on debug mode", default=False, required=False)
    return parser.parse_args(args)


def main(args):
    parsed_args = parse_args(args)

    if parsed_args.jobID is None:
        validator = Validator(debug_mode=parsed_args.d)
    else:
        validator = Validator(parsed_args.jobID, debug_mode=parsed_args.d)

    json_input_file = read_json_file(parsed_args.i, validator.project_logger,
                                     json_not_found_error_code="E.val.JSN.1",
                                     json_bad_format_error_code="E.val.JSN.2"
                                     )

    # The config file as a dictionary, with long key names replaced with their short names
    config_dict = validator.trim_config_file_keys(json_input_file)

    # The key type dictionary
    key_type_dict = read_json_file(parsed_args.KeyTypeFile, validator.project_logger,
                                   json_not_found_error_code="E.val.JSN.1",
                                   json_bad_format_error_code="E.val.JSN.2"
                                   )

    # Check the values present in the configuration file
    validator.validate_keys(config_dict, key_type_dict)


if __name__ == "__main__":
    main(sys.argv[1:])
