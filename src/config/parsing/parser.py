'''
    This script parses path information of Bioinformatics tools to a json Input File

                              Script Options

    -I      "The input Tool info text file "                                  (Required)
    -O      "json Input file to which path information is written into"       (Required)
'''

import ast
import re
import json
import argparse
import sys
import string


'''
 This function creates user-friendly command line interfaces. It defines what are the requires arguments for the script
'''
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        nargs=2,
                        metavar=('Tool_info_File', 'Inputs File'),
                        action='append',
                        help='Input Tools File which has the paths to all the executables',
                        required=True
                        )
    parser.add_argument("-o",
                        help='Output json File which has to be loaded with the path to executables',
                        required=True
                        )
    return parser.parse_args()


'''
 Comment line from the script are removed in this function
'''
def remove_comments(input_lines):
    filtered_lines = []
    for line in input_lines:
        if line != "":
            # If the first non-space character is a '#', exclude it
            if line.strip()[0] != '#':
                filtered_lines.append(line)
    return filtered_lines


'''
 This function definition is to create a Key:Value pair of Tools and Paths
'''
def create_key_value_pairs(tool_input_lines):

    key_value_pairs = []

    for line in tool_input_lines:
        if "=" not in line:
            sys.exit("ERROR: no equals sign present in the following line: " + line)
        else:
            # Get the position of the first equals sign
            split_pos = line.index("=")
            key = line[0:split_pos]
            value = line[split_pos + 1:]
            key_value_pairs.append([key, value])
    return key_value_pairs


'''
 In this function multiple checks are performed to see if the Tool:Path pairs are unique and void of typos
 
 The input is a List of Lists. The function returns separate lists for Tools and Paths'''
def create_tools_path_list(input_tools_path):

    # List comprehensions to create individual list for Tools and Paths
    Tools = list(map(lambda x: x[0], input_tools_path))
    Paths = list(map(lambda x: x[1], input_tools_path))

    for tool, executable in list(zip(Tools, Paths)):

        # Verifies if all the Tools have a corresponding Path to it
        if executable == '':
            sys.exit("ERROR: Tool " + tool + " is missing the path to the executable")

        # Check to see if the Paths are enclosed in double quotes in the Tool info text file
        elif executable[0] != '"' or executable[-1] != '"':
            sys.exit("ERROR: Executable for " + tool + " is not enclosed in double quotes")

        # Check to verify if special characters are present in the Paths information 
        for specialChar in string.punctuation.replace('/', '').replace('"', '').replace('-', '').replace('.', ''):
            if specialChar in executable:
                sys.exit("ERROR: Executable for " + tool + ' has an invalid special character: "' + specialChar + '"')
        
        # Check whether any key is present multiple times
        if Tools.count(tool) > 1:
            sys.exit("ERROR: Tool " + tool + " is listed twice in the tools list text file")

        else:
            return (Tools, Paths)

'''
 In this function Tools in the json input file is parsed with the Path information 
'''
def capture_executables(Output_Dict, Tools, Paths):
    for key, value in Output_Dict.items():
        for i in range(len(Tools)):
            # Regex check to find the list of tools in the json input file 
            if re.findall(Tools[i], key):

                # Removing double quotes from the paths and parsing path information to the json input file
                Paths[i] = Paths[i].replace('"', '')
                Output_Dict[key] = Paths[i]
    return Output_Dict


def main():

    # The argParser Function call
    args = parse_args()

    # List to hold all the input files
    list_of_input_files = []

    # This loop converts the list of lists into a list
    for InputFiles in args.i:
        for file in InputFiles:
            list_of_input_files.append(file)

    # This loops iterates through all the input files and parses the path information
    for i in range(len(list_of_input_files)):

        f1 = open(list_of_input_files[i], "r")
    
        # The lines in the file are split using line breaks
        ToolsTextLines = f1.read().splitlines()  
    
        # Filter out empty strings
        InputLines = list(filter(None, ToolsTextLines))

        # Function call remove_comments function
        InputLines = remove_comments(InputLines)
    
        # The create_key_value_pairs function call
        Tools_Paths_Input = create_key_value_pairs(InputLines)
        f1.close()
 
        # Function call to tools_and_paths function
        (Tools, Paths) = create_tools_path_list(Tools_Paths_Input)

        with open(args.o, "r") as OutputFile:
            Strings = OutputFile.read()

        OutputFile.close()

        # The "ast" python module helps convert a json file to a Python Dictionary
        OutputDict = ast.literal_eval(Strings)
 
        # The executableCapture function call
        OutputDict = capture_executables(OutputDict, Tools, Paths)

        # The "json" python module can be used to convert a Python Dictionary to a json file
        with open(args.o, "r+") as updated_json:
            json.dump(OutputDict, updated_json, indent=4)

        updated_json.close()


if __name__ == '__main__':
    main()

