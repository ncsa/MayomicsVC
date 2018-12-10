#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## deliver_alignment.sh MANIFEST, USAGE DOCS, SET CHECKS
#-------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST

*****************************************************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*****************************************************************************

MANIFEST
echo -e "${MANIFEST}"





read -r -d '' DOCS << DOCS

#############################################################################
#
# Deliver results of Design Block 1: trim-seq, align, sort, dedup. 
# Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 deliver_alignment.sh     -s           <sample_name>
                          -b           <aligned.sorted.deduped.bam>
                          -j           <WorkflowJSONfile>
                          -f           </path/to/delivery_folder>
                          -F           </path/to/shared_functions.sh>
                          -d           turn on debug mode

 EXAMPLES:
 deliver_alignment.sh -h     # get help message
 deliver_alignment.sh -s sample_name -b aligned.sorted.deduped.bam -j Workflow.json -f /path/to/delivery_folder -F /path/to/shared_functions.sh -d

#############################################################################

DOCS





set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=deliver_alignment.sh
SGE_JOB_ID=TBD   # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID






#-------------------------------------------------------------------------------------------------------------------------------
## LOGGING FUNCTIONS
#-------------------------------------------------------------------------------------------------------------------------------

function checkArg()
{
    if [[ "${OPTARG}" == -* ]]; then
        echo -e "\nError with option -${OPT} in command. Option passed incorrectly or without argument.\n"
        echo -e "\n${DOCS}\n"
        exit 1;
    fi
}






#-------------------------------------------------------------------------------------------------------------------------------
## GETOPTS ARGUMENT PARSER
#-------------------------------------------------------------------------------------------------------------------------------

## Check if no arguments were passed
if (($# == 0))
then
        echo -e "\nNo arguments passed.\n\n${DOCS}\n"
        exit 1
fi

## Input and Output parameters
while getopts ":hs:b:j:f:F:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo -e "\n${DOCS}\n"
			exit 0
                        ;;
                s )  # Sample name
                        SAMPLE=${OPTARG}
                        checkArg
                        ;;
                b )  # Full path to the input BAM file
                        BAM=${OPTARG}
			checkArg
                        ;;
                j )  # Full path to the workflow JSON file
                        JSON=${OPTARG}
                        checkArg
                        ;;
                f )  # Path to delivery folder
                        DELIVERY_FOLDER=${OPTARG}
                        checkArg
                        ;;
                F )  # Path to shared_functions.sh
                        SHARED_FUNCTIONS=${OPTARG}
                        checkArg
                        ;;
                d )  # Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d
                        echo -e "\nDebug mode is ON.\n"
			set -x
                        ;;
		\? )  # Check for unsupported flag, print usage and exit.
                        echo -e "\nInvalid option: -${OPTARG}\n\n${DOCS}\n"
                        exit 1
                        ;;
                : )  # Check for missing arguments, print usage and exit.
                        echo -e "\nOption -${OPTARG} requires an argument.\n\n${DOCS}\n"
                        exit 1
                        ;;
        esac
done





#-------------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS
#-------------------------------------------------------------------------------------------------------------------------------


source ${SHARED_FUNCTIONS}

## Check if Sample Name variable exists
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.deliver_alignment.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.deliver_alignment.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${BAM+x}" "Missing BAM option: -b" $LINENO
checkFile ${BAM} "Input BAM file ${BAM} is empty or does not exist." $LINENO
checkFile ${BAM}.bai "Input BAM index file ${BAM}.bai is empty or does not exist." $LINENO

checkVar "${JSON+x}" "Missing JSON option: -j" $LINENO
checkFile ${JSON} "Input JSON file ${JSON} is empty or does not exist." $LINENO

checkVar "${DELIVERY_FOLDER+x}" "Missing delivery folder option: -f" $LINENO





#-------------------------------------------------------------------------------------------------------------------------------
## MAKE DELIVERY FOLDER
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Creating the Delivery folder."

## Make delivery folder
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}: Creating Design Block 1 delivery folder. " ' INT TERM EXIT
makeDir ${DELIVERY_FOLDER} "Delivery folder ${DELIVERY_FOLDER}" ${TRAP_LINE}
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Created the Design Block 1 delivery folder."







#-------------------------------------------------------------------------------------------------------------------------------
## DELIVERY
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Copying Design Block 1 outputs into Delivery folder."

## Copy the files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}: Copying BAM into delivery folder. " ' INT TERM EXIT
cp ${BAM} ${DELIVERY_FOLDER}/${SAMPLE}.bam
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Aligned sorted deduped BAM delivered."


TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}: Copying BAM.BAI into delivery folder. " ' INT TERM EXIT
cp ${BAM}.bai ${DELIVERY_FOLDER}/${SAMPLE}.bam.bai
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Aligned sorted deduped BAM index delivered."

## Copy the JSON over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}: Copying JSON into delivery folder. " ' INT TERM EXIT
cp ${JSON} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Workflow JSON delivered."




#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index, and JSON. Open read permissions to the user group
checkFile ${DELIVERY_FOLDER}/${SAMPLE}.bam "Delivered BAM file ${DELIVERY_FOLDER}/${SAMPLE}.bam is empty" $LINENO
checkFile ${DELIVERY_FOLDER}/${SAMPLE}.bam.bai "Delivered BAM index file ${DELIVERY_FOLDER}/${SAMPLE}.bam.bai is empty" $LINENO

JSON_FILENAME=`basename ${JSON}` 
checkFile ${DELIVERY_FOLDER}/${JSON_FILENAME} "Delivered workflow JSON file ${DELIVERY_FOLDER}/${JSON_FILENAME} is empty" $LINENO

chmod g+r ${DELIVERY_FOLDER}/${SAMPLE}.bam
chmod g+r ${DELIVERY_FOLDER}/${SAMPLE}.bam.bai
chmod g+r ${DELIVERY_FOLDER}/${JSON_FILENAME}

logInfo "[DELIVERY] Alignment block delivered. Have a nice day."


## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
