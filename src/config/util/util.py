#!/usr/bin/env python
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("The script (" + sys.argv[0] + ") requires Python version 3.6 or higher")
    sys.exit(1)


'''
Takes in a list of input lines and returns that list of lines with comment lines removed

  Comment lines are any lines where the first non-whitespace character is a '#'
'''
def remove_comments(input_lines):
    def is_comment_line(line):
        # If the first non-space character is a '#', the line is a comment
        return True if line.strip()[0] == '#' else False
    #  This filters out all comment lines
    return list(filter(lambda x: not is_comment_line(x), input_lines))


'''
Takes a list of input lines and returns that list with any empty lines removed
'''
def remove_blank_lines(input_lines):
    #  This filters out all blank lines
    return list(filter(lambda x: x == "", input_lines))
