#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## deliver_haplotyperVC.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Deliver results of HaplotyperVC block: vcf for snps and indels, their index files, and workflow JSON. 
# Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 deliver_haplotyperVC.sh  -s           <sample_name>
                          -r           <OutputGVCF> 
                          -j           <WorkflowJSONfile>
                          -f           </path/to/delivery_folder>
                          -F           </path/to/shared_functions.sh>
                          -d           turn on debug mode

 EXAMPLES:
 deliver_haplotyperVC.sh -h
 deliver_haplotyperVC.sh -s sample_name -r Output.g.vcf -j Workflow.json -f /path/to/delivery_folder -F /path/to/shared_functions.sh -d

#############################################################################

DOCS





set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=deliver_haplotyperVC.sh
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
                r )  # Full path to the GVCF file
                        GVCF=${OPTARG}
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
ERRLOG=${SAMPLE}.deliver_haplotyperVC.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.deliver_haplotyperVC.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar ${GVCF} "Missing GVCF option: -r" $LINENO
checkFile ${GVCF} "Input GVCF file ${GVCF} is empty or does not exist" $LINENO
checkFile ${GVCF}.idx "Input GVCF index file ${GVCF}.idx is empty or does not exist" $LINENO

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
trap 'logError " $0 stopped at line ${TRAP_LINE}. Creating HaplotyperVC block delivery folder. " ' INT TERM EXIT
makeDir ${DELIVERY_FOLDER} "Delivery folder ${DELIVERY_FOLDER}" $LINENO
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[DELIVERY] Created the HaplotyperVC block delivery folder."







#-------------------------------------------------------------------------------------------------------------------------------
## DELIVERY
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Copying HaplotyperVC block outputs into Delivery folder."

## Copy the snp files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying GVCF into delivery folder. " ' INT TERM EXIT
cp ${GVCF} ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[DELIVERY] GVCF delivered."


TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying GVCF.IDX into delivery folder. " ' INT TERM EXIT
cp ${GVCF}.idx ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf.idx
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[DELIVERY] GVCF.IDX delivered."


## Copy the JSON over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying JSON into delivery folder. " ' INT TERM EXIT
cp ${JSON} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[DELIVERY] Workflow JSON delivered."



#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output GVCF and index, and JSON. Open read permissions to the user group
checkFile ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf "Delivered the GVCF file ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf is empty." $LINENO
checkFile ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf.idx "Delivered the GVCF index file ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf.idx is empty." $LINENO

JSON_FILENAME=`basename ${JSON}`
checkFile ${DELIVERY_FOLDER}/${JSON_FILENAME} "Delivered workflow JSON file ${DELIVERY_FOLDER}/${JSON_FILENAME} is empty" $LINENO

chmod g+r ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf
chmod g+r ${DELIVERY_FOLDER}/${SAMPLE}.g.vcf.idx
chmod g+r ${DELIVERY_FOLDER}/${JSON_FILENAME}

logInfo "[DELIVERY] HaplotyperVC block delivered. Have a nice day."


## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
