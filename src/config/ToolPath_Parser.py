# Import all the required modules
import ast
import re
import json
import argparse
import sys
import string

# Function definition 
def argParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input Tools File which has the paths to all the executables", required=True)
    parser.add_argument("-o", help="Output json File which has to be loaded with the path to executables", required=True)
    return parser.parse_args()


def commentRemoval(input_lines):                                                                                       
    filtered_lines = []                                                                                                 
    for line in input_lines:                                                                                            
        if line != "":                                                                                               
           # If the first non-space character is a '#', exclude it
           if line.strip()[0] != '#':                                                                                   
               filtered_lines.append(line)
    return filtered_lines


def keyValuePairCreation( toolInputLines ):
    keyValuePairs = []
    for line in toolInputLines:
        if "=" not in line:
            sys.exit("ERROR: no equals sign present in the following line: " + line)
        else:
            # Get the position of the first equals sign
            split_pos = line.index("=")
            key = line[0:split_pos]
            value = line[split_pos + 1:]
            keyValuePairs.append([key, value])
    return keyValuePairs


def toolsandPaths( Input_Tools_Path ):
    Tools = list(map(lambda x: x[0], Input_Tools_Path))
    Paths = list(map(lambda x: x[1], Input_Tools_Path))

    for tool, executable in list(zip(Tools, Paths)):
        if executable == '':
            sys.exit("ERROR: Tool " + tool + " is missing the path to the executable")

        elif executable[0] != '"' or executable[-1] != '"':
            sys.exit("ERROR: Executable for " + tool + " is not enclosed in double quotes")

        for specialChar in string.punctuation.replace('/', '').replace('"', '').replace('-','').replace('.',''):
            if specialChar in executable:
                sys.exit("ERROR: Executable for " + tool + ' has an invalid special character: "' + specialChar + '"')
        
        # Check whether any key is present multiple times
        if (Tools.count(tool) > 1):
            sys.exit("ERROR: Tool " + tool + " is listed twice in the tools list text file")

        else:
            return (Tools, Paths)


def executablesCapture( Output_Dict, Tools, Paths ):

    for key, value in Output_Dict.items():

        for i in range(len(Tools)):
            if (re.findall(Tools[i], key)):
                Paths[i] = Paths[i].replace('"','')
                Output_Dict[key] = Paths[i]

    return Output_Dict


def main():   
    args = argParser() 
    
    InputFile = open(args.i,"r")
    
    ToolsTextLines = InputFile.read().splitlines()  
    
    # Filter out empty strings
    InputLines = list(filter(None, ToolsTextLines))

    InputLines = commentRemoval(InputLines)
    
    Tools_Paths_Input = keyValuePairCreation(InputLines)
    InputFile.close()
 
    (Tools, Paths) = toolsandPaths(Tools_Paths_Input)

    with open(args.o,"r") as OutputFile:
        Strings = OutputFile.read()

    OutputFile.close()

    OutputDict = ast.literal_eval(Strings)

    OutputDict = executablesCapture(OutputDict, Tools, Paths)

    with open(args.o, "r+") as updated_json:
        json.dump(OutputDict, updated_json, indent=4)

    updated_json.close()

if __name__ == '__main__':
    main()

