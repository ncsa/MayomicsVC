import ast
import re
import json
import argparse

'''                                                                                                                     
JRH comments:                                                                                                           
                                                                                                                        
The "main" function should be called main                                                                               
                                                                                                                        
'''                                                                                                                     
                                                                                                                        
                                                                                                                        
                                                                                                                        
def ArgParser():                                                                                                        
    parser = argparse.ArgumentParser()                                                                                  
    parser.add_argument("-i", help="Input Tools File which has the paths to all the executables", required=True)        
    parser.add_argument("-o", help="Output json File which has to be loaded with the path to executables", required=True)
    return parser.parse_args()                                                                                          
                                                                                                                        
                                                                                                                        
                                                                                                                        
'''                                                                                                                     
                                                                                                                        
                                                                                                                        
'''                                                                                                                     
                                                                                                                        
## fthjsryi                                                                                                             
        #dghj fghj f                                                                                                    
                                                                                                                        
'''                                                                                                                     
Description of func here                                                                                                
                                                                                                                        
'''                                                                                                                     
def comment_removal(input_lines):                                                                                       
    filtered_lines = []                                                                                                 
    for line in input_lines:                                                                                            
        if line != "":                                                                                                  
           # If the first non-space character is a '#', exclude it                                                      
           if line.strip()[0] != '#':                                                                                   
               filtered_lines.append(line)                                                                              
    return filtered_lines

def CommentsRemoval ( InputLines ):
    for ListElement in InputLines:
        for letter in range(len(ListElement)):
            if (ListElement[letter] == '#'):
                InputLines.remove(ListElement)
    return InputLines

'''
What happens if a line has no '=' sign?

Currently: It will not split, and the list of Key-value pairs will have an element with only 1 item, which breaks
      the code later

We want an intelligent error message

Option 1: Ignore the line, maybe throw the user a message printing the line and saying it was ignore

Option 2: Throw a fatal error, printing the line and explaining

Lets go with option 2

What if there is more than 1 equals sign (say, a path has an equals sign, or some string contains one)?

Currently:
    The key value list that is included as an element of the return list will have more than 2 entries. Since
    the later code only looks at the first two elements, this will not break the parser. But the value that the user
    expected to have and what they will get will be different.

    For example:    lets say the line was 'KEY="DUMP=VALUE"', the current code will see the following key and values

'''



def KeyValuePairCreation( ToolsInputLines ):
    Input_Tools_Paths = []
    for i in ToolsInputLines:
        KV_pair = i.split("=")
        Input_Tools_Paths.append(KV_pair)
    return Input_Tools_Paths

def ToolsandPaths( Input_Tools_Path ):

    Tools = list(map(lambda x: x[0], Input_Tools_Path))
    Paths = list(map(lambda x: x[1], Input_Tools_Path))
    return (Tools, Paths)

def ExecutablesCapture( Output_Dict, Tools, Paths ):

    for key, value in Output_Dict.items():

        for i in range(len(Tools)):
            if (re.findall(Tools[i], key)):
                Paths[i] = Paths[i].replace('"','')
                Output_Dict[key] = Paths[i]

    return Output_Dict

def Tool_Parser():   
    args = ArgParser() 
    InputFile = open(args.i,"r")
    ToolsTextLines = InputFile.read().splitlines()
    
    InputLines = list(filter(None, ToolsTextLines))
 
    InputLines = CommentsRemoval(InputLines)

    Tools_Paths_Input = KeyValuePairCreation(InputLines)
    InputFile.close()
 
    (Tools, Paths) = ToolsandPaths(Tools_Paths_Input)

    with open(args.o,"r") as OutputFile:
        Strings = OutputFile.read()

    OutputFile.close()

    OutputDict = ast.literal_eval(Strings)

    OutputDict = ExecutablesCapture(OutputDict, Tools, Paths)

    with open(args.o, "r+") as updated_json:
        json.dump(OutputDict, updated_json, indent=4)

    updated_json.close()

if __name__ == '__main__':
    Tool_Parser()

