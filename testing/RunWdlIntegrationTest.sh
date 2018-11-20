#!/bin/bash

# EXAMPLE USAGE 
# ./AutomatingTesting.sh /path/to/test/folder src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl "-i Config/memory_info.txt -i Config/run_info.txt -i Config/sample_info.txt -i Config/tool_info.txt"



set -x
set -o errexit
set -o pipefail
set -o nounset

# ASSUMPTIONS
# cd into your test working folder
# in there do git pull on the repo, which creates the MayomicsVC/ folder
# also in there have folders Config/, Jsons/, Inputs/ 
# - although it is also a good practice to keep all input data in one central place, not in test folder
# Let's name that workspace ${TestingMayomics}
# TestingMayomics="/projects/mgc/Project_1/LSM/CleaningUpWDLtails/Multilane/"
# 

#Define the path to test folder
TestingMayomics=$1

#define which stage of the workflow is being tested: relative path to .wdl file within the MayomicsVC repository
WorkflowBeingTested=$2;
BaseNameOfWorkflowBeingTested=`basename ${WorkflowBeingTested} .wdl`

#define which config files to use in this test
#should be in format "-i /full/path/to/config1.txt -i /full/path/to/config2.txt etc..."
#or can automatically grab all configs in a folder, if you like - easy enough to code that up
ConfigsBeingUsed=$3;



#HERE BEGINS THE ACTUAL CODE THAT CAN BE RUNNABLE FROM THIS SCRIPT OR COPY-PASTED

cd ${TestingMayomics};

# these are specific to our system; should we perhaps make these paths be input variables?
source /etc/profile.d/modules.sh
module load cromwell/cromwell-34;
module load python/python-3.6.1;



#### JSON STUFF ##########################
# create JSON template
cd "./MayomicsVC/"; 
java -jar ${WOMTOOL} inputs ${WorkflowBeingTested} > ../Jsons/${BaseNameOfWorkflowBeingTested}.json;
cd ../;

#populate the JSON template
python MayomicsVC/src/config/config_parser.py ${ConfigsBeingUsed} --jsonTemplate Jsons/${BaseNameOfWorkflowBeingTested}.json -o Jsons/${BaseNameOfWorkflowBeingTested}.FilledIn.json;

#validate the JSON template
python MayomicsVC/src/config/key_validation.py -i Jsons/${BaseNameOfWorkflowBeingTested}.FilledIn.json --KeyTypeFile MayomicsVC/src/config/key_types.json;




#####    create the zip file ##############
cd MayomicsVC ; 
zip -r MayomicsVC.zip ./ ;
mv MayomicsVC.zip ../ ;
cd ../ ;




##### RUN THE WORKFLOW #############

java -jar ${CROMWELL} run ./MayomicsVC/${WorkflowBeingTested} -i Jsons/${BaseNameOfWorkflowBeingTested}.FilledIn.json -p MayomicsVC.zip ;




