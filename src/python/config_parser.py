#!/usr/bin/env python3


import argparse
import sys
from config.parser.parsing import Parser


def parse_args(args):
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
    return parser.parse_args(args)


def main(args):
    parsed_args = parse_args(args)

    # Instantiation of the Parser class
    if parsed_args.jobID is None:
        k_v_parser = Parser(debug_mode=parsed_args.d)
    else:
        k_v_parser = Parser(parsed_args.jobID, debug_mode=parsed_args.d)

    # Fill in the json template file with values from the Key="Value" formatted input files
    k_v_parser.fill_in_json_template(parsed_args.i, parsed_args.jsonTemplate, parsed_args.o)


if __name__ == '__main__':
    main(sys.argv[1:])
