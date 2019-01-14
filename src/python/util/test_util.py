#!/usr/bin/env python

import unittest
import os
from util.log import ProjectLogger
from util.util import read_json_file
from util.util import flatten

# This is the full path to the validation package
package_full_path = os.path.abspath(os.path.dirname(__file__))


class TestUtil(unittest.TestCase):
    # Create a project logger just for these unit tests
    project_logger = ProjectLogger(name="util-test", job_id="NA")

    # Turn the project logger off during UnitTesting, so the end user is not confused by error messages
    #   (Some tests are designed to fail, so they will log "ERROR" messages that are expected)
    project_logger.logger.disabled = True

    def test_read_json_file(self):
        expected_key_types_dict = {"key_1": "Value_1", "key_2": "Value_2", "key_3": "Value_3"}

        try:
            json_dict = read_json_file(
                package_full_path + "/test_resources/test.json",
                project_logger=self.project_logger,
                json_bad_format_error_code="BadFormat",
                json_not_found_error_code="NotFound"
            )
            self.assertEqual(json_dict, expected_key_types_dict)
        # This catch safely handles the function failing and prevents python itself from exiting (The function will
        #   fail if the input json cannot be found or if the file is improperly formatted)
        except SystemExit:
            self.fail(msg="'read_json_file' tried to exit. The unit test 'test_read_json_file' must fail")

    '''
    This method tries to find a file that does not exist; it is supposed to fail
    '''
    def test_read_json_file_path_failure(self):
        with self.assertRaises(SystemExit):
            read_json_file("/this/path/does/not/exist.json", project_logger=self.project_logger,
                           json_bad_format_error_code="BadFormat", json_not_found_error_code="NotFound"
                           )

    def test_flatten(self):
        input_list = [[1,2], [9,[[11, 4], 87, "a"]]]
        expected = [1, 2, 9, 11, 4, 87, "a"]
        actual = flatten(input_list)

        self.assertEqual(expected, actual)

