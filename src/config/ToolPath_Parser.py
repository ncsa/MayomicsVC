import ast
import re
import json

with open("/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/config/Test_Tools.txt","r") as InputFile:
    InputLines = InputFile.read().splitlines()

Tools_Paths_Input = []

for e in InputLines:
    KV_pair = e.split("=")
    Tools_Paths_Input.append(KV_pair)

InputFile.close()

Tools = list(map(lambda x: x[0], Tools_Paths_Input))
Paths = list(map(lambda x: x[1], Tools_Paths_Input))


with open("/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/config/Output_File.json","r") as OutputFile:
    String = OutputFile.read()

Tools_Path_Output = []

OutputFile.close()

OutputDict = ast.literal_eval(String)


for key, value in OutputDict.items():

    for i in range(len(Tools)):
        if (re.findall(Tools[i], key)):
            OutputDict[key] = Paths[i]

with open("/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/config/Output_File.json", "r+") as updated_json:
    json.dump(OutputDict, updated_json, indent=4)

updated_json.close()



