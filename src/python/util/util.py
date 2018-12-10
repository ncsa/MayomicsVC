#!/usr/bin/env python

import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("The script (" + sys.argv[0] + ") requires Python version 3.6 or higher")
    sys.exit(1)

import json

def flatten(input_container):
    """
    :param container: collection of lists and tuples with arbitrary nesting levels
    :return: A 1-D list of the input's contents
    """
    def __flattener(container):
        for i in container:
            if isinstance(i, (list, tuple)):
                for j in __flattener(i):
                    yield j
            else:
                yield i

    return list(__flattener(input_container))


def read_json_file(json_file, project_logger, json_not_found_error_code, json_bad_format_error_code):
    """
    Reads a JSON formatted file, and returns its contents as a Python dictionary

    Inputs:
        The json file
        The project logger (for error messages to get passed to)
        The error code if a JSON file is not found
        The error code if a JSON file is improperly formatted
    """
    try:
        # Open the file
        with open(json_file) as F:
            # Read the file's contents as a string
            json_str = F.read()
            # Return the data as a Python dictionary
            return json.loads(json_str)
    except FileNotFoundError:
        project_logger.log_error(
            json_not_found_error_code,
            'Could not open json file "' + str(json_file) + '": JSON file could not be found'
        )
        sys.exit(1)
    except json.decoder.JSONDecodeError:
        project_logger.log_error(
            json_bad_format_error_code,
            'Could not open json file "' + str(json_file) + '": JSON file not formatted properly'
        )
        sys.exit(1)
