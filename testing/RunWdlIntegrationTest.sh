#!/bin/bash

# EXAMPLE USAGE 
# ./AutomatingTesting.sh /projects/mgc/Project_1/LSM/CleaningUpWDLtails/Multilane src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl "-i Config/memory_info.txt -i Config/run_info.txt -i Config/sample_info.txt -i Config/tool_info.txt"



set -x
set -o errexit
set -o pipefail
set -o nounset

# ASSUMPTIONS
# cd into your test working folder
# in there do git pull on the repo, which creates the MayomicsVC/ folder
# also in there have folders Config/, Jsons/, Inputs/ - although it is also a good practice  to keep all input data in one central place, not in test folder
# Let's name that workspace ${TestingMayomics}
# TestingMayomics="/projects/mgc/Project_1/LSM/CleaningUpWDLtails/Multilane/"
# Folder that includes hg38 simurlated data: /projects/bioinformatics/DataPacks/human/Hg38_chr20_21_22_simulated_data_Jan_2018/fastqs/
# 

#Define the path to test folder
TestingMayomics=$1

#define which stage of the workflow is being tested: relative path to .wdl file within the MAyomicsVC repository
WorkflowBeingTested=$2;
BaseNameOfWorkflowBeingTested=`basename ${WorkflowBeingTested} .wdl`

#define which config files to use in this test
#should be in format "-i /full/path/to/config1.txt -i /full/path/to/config2.txt -i /full/path/to/config3.txt etc..."
#or can automatically grab all configs in a folder, if you like - easy enough to code that up
ConfigsBeingUsed=$3;



#HERE BEGINS THE ACTUAL CODE THAT CAN BE RUNNABLE FROM THIS SCRIPT OR COPY-PASTED

cd ${TestingMayomics};

# these are specific to our system and I'd rather they were not explicitly mentioned in a script that hangs on github
# should we perhaps make these paths be input variables?
# module load directly does not work from within Bash, and the workaround seems complicated; I defined the variables directly for that reason
source /etc/profile.d/modules.sh
module load /usr/local/apps/bioapps/modules/cromwell/cromwell-34;
module load python/python-3.6.1;
#CROMWELL=/usr/local/apps/bioapps/cromwell/cromwell-34/cromwell-34.jar
#WOMTOOL=/usr/local/apps/bioapps/cromwell/cromwell-34/womtool-34.jar
#export PATH=/usr/local/apps/bioapps/python/Python-3.6.1:/usr/local/apps/bioapps/python/Python-3.6.1/bin:${PATH}
#export LD_LIBRARY_PATH=/usr/local/apps/bioapps/python/Python-3.6.1/lib:${PATH}



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




