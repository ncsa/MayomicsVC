#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## deliver_somaticVC.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Deliver results of SomaticVC block: vcf for snps and indels, their index files, and workflow JSON. 
# Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 deliver_somaticVC.sh     -s           <sample_name>
                          -r           <ResultantVcf> 
                          -j           <WorkflowJSONfile>
                          -f           </path/to/delivery_folder>
                          -F           </path/to/shared_functions.sh>
                          -d           turn on debug mode

 EXAMPLES:
 deliver_somaticVC.sh -h
 deliver_somaticVC.sh -s sample_name -r ResultantVcf.vcf -j Workflow.json -f /path/to/delivery_folder -F /path/to/shared_functions.sh -d

#############################################################################

DOCS





set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=deliver_somaticVC.sh
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
while getopts ":hs:r:j:f:F:d" OPT
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
                r )  # Full path to the VCF file
                        VCF=${OPTARG}
                        checkArg
                        ;;
                j )  # Full path to the workflow JSON file
                        JSON=${OPTARG}
                        checkArg
                        ;;
                f)   # Path to delivery folder
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
checkVar ${SAMPLE} "Missing sample name option: -s" $LINENO

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.deliver_somaticVC.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.deliver_somaticVC.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar ${VCF} "Missing VCF option: -r" $LINENO
checkFile ${VCF} "Input VCF file ${VCF} is empty or does not exist" $LINENO
checkFile ${VCF}.idx "Input VCF index file ${VCF}.idx is empty or does not exist" $LINENO

checkVar ${JSON} "Missing JSON option: -j" $LINENO
checkFile ${JSON} "Input JSON file ${JSON} is empty or does not exist." $LINENO

checkVar ${DELIVERY_FOLDER} "Missing delivery folder option: -f" $LINENO





#-------------------------------------------------------------------------------------------------------------------------------
## MAKE DELIVERY FOLDER
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Creating the Delivery folder."

## Make delivery folder
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Creating SomaticVC block delivery folder. " ' INT TERM EXIT
makeDir ${DELIVERY_FOLDER} "Delivery folder ${DELIVERY_FOLDER}" $LINENO
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Created the SomaticVC block delivery folder."







#-------------------------------------------------------------------------------------------------------------------------------
## DELIVERY
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Copying SomaticVC block outputs into Delivery folder."

## Copy the snp files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying VCF into delivery folder. " ' INT TERM EXIT
cp ${VCF} ${DELIVERY_FOLDER}/${SAMPLE}.vcf
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Recalibrated VCF delivered."


TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying VCF.IDX into delivery folder. " ' INT TERM EXIT
cp ${VCF}.idx ${DELIVERY_FOLDER}/${SAMPLE}.vcf.idx
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Recalibrated VCF.IDX delivered."


## Copy the JSON over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying JSON into delivery folder. " ' INT TERM EXIT
cp ${JSON} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE}
logInfo "[DELIVERY] Workflow JSON delivered."



#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output VCF and index, and JSON. Open read permissions to the user group
checkFile ${DELIVERY_FOLDER}/${SAMPLE}.vcf "Delivered recalibrated VCF file ${DELIVERY_FOLDER}/${SAMPLE}.vcf is empty." $LINENO
checkFile ${DELIVERY_FOLDER}/${SAMPLE}.vcf.idx "Delivered recalibrated VCF index file ${DELIVERY_FOLDER}/${SAMPLE}.vcf.idx is empty." $LINENO

JSON_FILENAME=`basename ${JSON}`
checkFile ${DELIVERY_FOLDER}/${JSON_FILENAME} "Delivered workflow JSON file ${DELIVERY_FOLDER}/${JSON_FILENAME} is empty" $LINENO

chmod g+r ${DELIVERY_FOLDER}/${SAMPLE}.vcf
chmod g+r ${DELIVERY_FOLDER}/${SAMPLE}.vcf.idx
chmod g+r ${DELIVERY_FOLDER}/${JSON_FILENAME}

logInfo "[DELIVERY] SomaticVC block delivered. Have a nice day."


## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
