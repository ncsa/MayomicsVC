#!/usr/bin/env python

import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("This script (" + sys.argv[0] + ") requires Python version 3.6 or higher")
    sys.exit(1)

import os
import argparse
import logging
import pathlib
from util.util import flatten
from util.log import ProjectLogger
from config.util.special_keys import OPTIONAL_KEYS

"""
Exit code Rules:

1. Exit codes in this module are only given when an error has occurred, so they will all start with 'E.'
2. The letters 'val.' because they are coming from the validation component of the code
3. A three letter code that hints at what the problem was
4. A number that can help to differentiate similar error codes


Error Code List
========================================================================================================================
E.val.JSN.1 = An input JSON file could not be found
E.val.JSN.2 = An input JSON file was not formatted properly

E.val.OuD.1 = An output directory already exists

E.val.ROD.1 = A directory that is supposed to be readable and executable could not be found
E.val.ROD.2 = A read-only directory does not have read permissions for the current user
E.val.ROD.3 = A read-only directory does not have executable permissions for the current user


E.val.ExF.1 = A file expected to be executable could not be found
E.val.ExF.2 = A file expected to be executable was not executable by the current user

E.val.Fil.1 = A regular file could not be found
E.val.Fil.2 = A regular file was not readable by the current user

E.val.Boo.1 = A value expected to be a boolean type was not
E.val.Int.1 = A value expected to be an integer was not
E.val.Dec.1 = A value expected to be a Decimal (float or double) was not
E.val.UNK.1 = A type listed in the key-types file was not recognized as a valid type

E.val.Bug.1 = The DebugMode key had a non-empty value that was not '-d'
"""


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


class Validator:
    def __init__(self, job_id="NA", debug_mode=False):
        # Initialize the project logger (with the level at INFO by default and DEBUG if in debug mode)
        if debug_mode:
            self.project_logger = ProjectLogger(job_id, "validator.key_validation.Validator", logging.DEBUG)
        else:
            self.project_logger = ProjectLogger(job_id, "validator.key_validation.Validator")
        self.job_id = job_id
        self.debug_mode = debug_mode

    @staticmethod
    def trim_config_file_keys(long_key_dict):
        """
        In the json config file, the keys are in a long format, such as MainTask.MinorTask.KeyName.
            We want the last item only: the KeyName

        This function trims off all but the last name (we can assume there are no '.' characters in the KeyName)
        The input is the key-value dict with the original long keys, and the output is the same dict with short keys
        """
        short_key_dict = {}
        for long_key in long_key_dict.keys():
            # Split the long key by '.' and select the last string (if no '.' present, the full string is selected)
            short_key = long_key.split('.')[-1]
            # New dictionary with new key that points to the original value
            short_key_dict[short_key] = long_key_dict[long_key]
        return short_key_dict

    @staticmethod
    def __file_exists(file_path):
        """
        Confirms that a file exists
        Returns a boolean
        """
        path = pathlib.Path(file_path)
        return path.is_file()

    def __file_is_readable(self, file_path):
        """
        Confirms that a file exists and is readable

        Returns one of three possible strings: (Success, FileNotFound, or FileNotReadable)
        """
        if self.__file_exists(file_path):
            # Check to see if the file on this path is readable by the user calling this script
            is_readable = os.access(file_path, os.R_OK)
            if is_readable:
                return "Success"
            else:
                return "FileNotReadable"
        else:
            return "FileNotFound"

    def __file_is_executable(self, file_path):
        """
        Confirms that a file exists and is executable

        Returns one of three possible strings: (Success, FileNotFound, or FileNotExecutable)
        """
        if self.__file_exists(file_path):
            # Check to see if the file on this path is executable by the user calling this script
            is_executable = os.access(file_path, os.X_OK)
            if is_executable:
                return "Success"
            else:
                return "FileNotExecutable"
        else:
            return "FileNotFound"

    @staticmethod
    def __directory_exists(directory_path):
        """
        Determines whether a given directory exists

        Returns a boolean
        """
        return os.path.isdir(directory_path)

    @staticmethod
    def __is_integer(input_string):
        """
        Will return true for integer strings like "9", "2", etc. but false for strings like "3.14" or "NotAString"
        """
        try:
            int(input_string)
            return True
        except ValueError:
            return False

    @staticmethod
    def __is_float(input_string):
        try:
            float(input_string)
            return True
        except ValueError:
            return False

    def check_key(self, key_name, key_value, key_type):
        """
        For a given key, confirm that its value is of the type that is expected; returns True if the key is
            valid and False if it is not

        If a value is validated, print a DEBUG message stating that X key was validated; return true
        If a value is of a type that has no validation defined (such as String), print a DEBUG message; return true
        If a faulty value is found, print an ERROR message; return false
        """
        lowered_key_type = key_type.lower()

        def make_message(message):
            return 'Input variable "' + key_name + '" points to "' + key_value + '", which ' + message

        # OPTIONAL_KEYS ###
        # If the key has an empty value, and it is in the OPTIONAL_KEYS list (see src/config/util/special_keys.py)
        #   then it is allowed to be blank, it counts as a valid entry
        if key_value == "" and key_name.lower() in OPTIONAL_KEYS:
            self.project_logger.log_debug(
                "The key '" + key_name + "' has an empty value; because it is an optional key, this is still valid"
            )
            # Stop the function here; do not try to validate this key
            return True

        # Directory ###
        if lowered_key_type in ("directory", "dir"):
            # Checks if the directory exists, and does not check permissions
            if not self.__directory_exists(key_value):
                self.project_logger.log_error("E.val.Dir.1", "The directory: '" + key_value + "' could not be found.")
                return False
            else:
                return True

        # Output Directory ###
        elif lowered_key_type in ("output_directory", "output_dir", "outputdir", "outputdirectory", "output directory"):
            # Only verifies that the directory does not already exist
            if self.__directory_exists(key_value):
                self.project_logger.log_error("E.val.OuD.1", "The output directory '" + key_value + "' already exists.")
                return False
            else:
                return True

        # Read-only Directory ###
        elif lowered_key_type in (
                "read-only-directory", "read_only_directory", "read-only_directory", "read-only",
                "readonly", "readonlydirectory", "readonlydir", "read-only-dir", "read_only_dir", "read only",
                "read only dir", "read only directory"
        ):
            # Checks if the directory exists, and has executable (ability to traverse into) and read (read contents).
            if not self.__directory_exists(key_value):
                self.project_logger.log_error("E.val.ROD.1", "The directory: '" + key_value + "' could not be found.")
                return False
            else:
                readable = os.access(key_value, os.R_OK)
                executable = os.access(key_value, os.X_OK)
                if not readable:
                    self.project_logger.log_error("E.val.ROD.2", "The directory: '" + key_value +
                                                  "' does not have read permissions for the current user."
                                                  )
                    return False
                elif not executable:
                    self.project_logger.log_error("E.val.ROD.3", "The directory: '" + key_value +
                                                  "' does not have executable permissions for the current user."
                                                  )
                    return False
                else:
                    # The directory was found and is readable and executable
                    return True

        # ExecFile ###
        elif lowered_key_type in ("execfile", "exec_file", "executable"):
            # exec_status can only be one of three strings: Success, FileNotFound, or FileNotExecutable
            exec_status = self.__file_is_executable(key_value)
            if exec_status == "Success":
                self.project_logger.log_debug(make_message('was found and is executable'))
                return True
            elif exec_status == "FileNotFound":
                self.project_logger.log_error("E.val.ExF.1", make_message('could not be found'))
                return False
            elif exec_status == "FileNotExecutable":
                self.project_logger.log_error("E.val.ExF.2", make_message('is not executable by the current user'))
                return False

        # File ###
        elif lowered_key_type == "file":
            # readable_status can only be one of three strings: Success, FileNotFound, or FileNotReadable
            readable_status = self.__file_is_readable(key_value)
            if readable_status == "Success":
                self.project_logger.log_debug(make_message('was found and is readable'))
                return True
            elif readable_status == "FileNotFound":
                self.project_logger.log_error("E.val.Fil.1", make_message('could not be found'))
                return False
            elif readable_status == "FileNotReadable":
                self.project_logger.log_error("E.val.Fil.2", make_message('is not readable by the current user'))
                return False

        # Boolean ###
        elif lowered_key_type in ("boolean", "bool"):
            # The key_value was not .lowered() because this value is case sensitive
            if key_value in ("true", "false"):
                self.project_logger.log_debug(make_message('is a valid boolean value'))
                return True
            else:
                self.project_logger.log_error(
                    "E.val.Boo.1",
                    make_message('is not a valid boolean value; enter either "true" or "false" (case sensitive)')
                )
                return False

        # String ###
        elif lowered_key_type in ("string", "str"):
            self.project_logger.log_debug(
                'Input variable "' + key_name + '" is a String type; no pre-workflow validation can take place'
            )
            return True

        # Integer ###
        elif lowered_key_type in ("integer", "int"):
            if self.__is_integer(key_value):
                self.project_logger.log_debug(make_message('is a valid integer'))
                return True
            else:
                self.project_logger.log_error('E.val.Int.1', make_message('is not a valid integer'))
                return False

        # Decimal ###
        elif lowered_key_type == "decimal":
            if self.__is_float(key_value):
                self.project_logger.log_debug(make_message('is a valid number'))
                return True
            else:
                self.project_logger.log_error('E.val.Dec.1', make_message('is not a valid number'))
                return False
        # DebugMode (special case where the only acceptable value is '-d')
        elif lowered_key_type == "debugmode":
            if key_value == "-d":
                self.project_logger.log_debug(make_message('is the correct debug flag'))
                return True
            else:
                self.project_logger.log_error(
                    'E.val.Bug.1', make_message("is not valid: DebugMode must be blank or '-d'")
                )
                return False
        # Other ###
        # (kill the program if an unknown type is provided in the type file)
        else:
            self.project_logger.log_error(
                'E.val.UNK.1',
                'Input variable "' + key_name + '" has the type "' + key_type +
                '" in the key types file, which is not a recognized type ' +
                '(see src/config/validation/README.md for a list of valid types)'
            )

            return False


    def validate_keys(self, configuration_dict, key_types_dict):
        """
        Loop through all keys and validate the types of all keys that have types listed in the Key-Types file
        """
        unchecked_keys = []
        checked_keys = []
        for key in configuration_dict:
            # Sort keys into groups of checked and unchecked keys
            checked_keys.append(key) if key in key_types_dict else unchecked_keys.append(key)

        # Log a warning message for all of the keys that will not be checked
        #   Sorted to ensure repeated problematic runs fail on the same key
        for key in sorted(unchecked_keys):
            self.project_logger.log_warning(
                'The variable "' + key + '" is present in the json config file but not in the key types file. ' +
                'Its value is unchecked! To check this value, add the key to the key type file'
            )
        # For all keys that have typing information, confirm that its value is valid
        #   Sorted to ensure repeated problematic runs fail on the same key
        for key in sorted(checked_keys):

            # If the type has square brackets, it is an array (or array of arrays, etc.)
            if key_types_dict[key][0] == "[" and key_types_dict[key][-1] == "]":
                trimmed_type = key_types_dict[key][1:-1]
                for entry in flatten(configuration_dict[key]):
                    key_is_valid = self.check_key(key, entry, trimmed_type)
                    if not key_is_valid:
                        sys.exit(1)
            # Check all other types of values
            else:
                key_is_valid = self.check_key(key, configuration_dict[key], key_types_dict[key])
                if not key_is_valid:
                    # The check key function itself wrote to the log; here we just exit because we have found an invalid key
                    sys.exit(1)

        # Write a success message to the log
        self.project_logger.log_info(
            'Key validation finished successfully with ' + str(self.project_logger.warnings_issued) +
            ' warning(s) issued'
        )


