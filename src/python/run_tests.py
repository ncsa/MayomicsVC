#!/usr/bin/env python3

import unittest

from config.parser.test_parsing import TestParser
from config.validator.test_key_validation import TestValidator
from util.test_util import TestUtil

"""
Each of these imports are Classes that inherit from unittest and contain unit tests (methods beginning with "test_")
  This script will automatically run all of the unit tests in those scripts 
"""

if __name__ == '__main__':
    unittest.main()
