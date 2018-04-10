import ast
import re
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input Tools File which has the paths to all the executables", required=True)
parser.add_argument("-o", help="Output json File which has to be loaded with the path to executables", required=True)
args =  parser.parse_args()

with open(args.i,"r") as InputFile:
    InputLines = InputFile.read().splitlines()

Tools_Paths_Input = []

for e in InputLines:
    KV_pair = e.split("=")
    Tools_Paths_Input.append(KV_pair)

InputFile.close()

Tools = list(map(lambda x: x[0], Tools_Paths_Input))
Paths = list(map(lambda x: x[1], Tools_Paths_Input))


with open(args.o,"r") as OutputFile:
    String = OutputFile.read()

Tools_Path_Output = []

OutputFile.close()

OutputDict = ast.literal_eval(String)


for key, value in OutputDict.items():

    for i in range(len(Tools)):
        if (re.findall(Tools[i], key)):
            Paths[i] = Paths[i].replace('"','')
            OutputDict[key] = Paths[i]

with open(args.o, "r+") as updated_json:
    json.dump(OutputDict, updated_json, indent=4)

updated_json.close()




if __name__ == "__main__":
    main()

