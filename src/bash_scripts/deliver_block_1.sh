#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## deliver_block_1.sh MANIFEST, USAGE DOCS, SET CHECKS
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
 deliver_block_1.sh       -b           <aligned.sorted.dedupped.bam>
                          -j           <WorkflowJSONfile>
                          -f           </path/to/delivery_folder>
                          -d           turn on debug mode

 EXAMPLES:
 deliver_block_1.sh -h     # get help message
 deliver_block_1.sh -b aligned.sorted.dedupped.bam -j Workflow.json -f /path/to/delivery_folder -d

#############################################################################

DOCS








set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=deliver_block_1.sh
SGE_JOB_ID=TBD  # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## LOGGING FUNCTIONS
#-------------------------------------------------------------------------------------------------------------------------------

## Logging functions
# Get date and time information
function getDate()
{
    echo "$(date +%Y-%m-%d'T'%H:%M:%S%z)"
}

# This is "private" function called by the other logging functions, don't call it directly,
# use logError, logWarn, etc.
function _logMsg () {
    echo -e "${1}"

    if [[ -n ${ERRLOG-x} ]]; then
        echo -e "${1}" | sed -r 's/\\n//'  >> "${ERRLOG}"
    fi
}

function logError()
{
    local LEVEL="ERROR"
    local CODE="-1"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    >&2 _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"

    if [[ -z ${EXITCODE+x} ]]; then
        EXITCODE=1
    fi

    exit ${EXITCODE};
}

function logWarn()
{
    local LEVEL="WARN"
    local CODE="0"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"
}

function logInfo()
{
    local LEVEL="INFO"
    local CODE="0"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"
}

function checkArg()
{
    if [[ "${OPTARG}" == -* ]]; then
        echo -e "\nError with option -${OPT} in command. Option passed incorrectly or without argument.\n"
        echo -e "\n${DOCS}\n"
        exit 1;
    fi
}

#-------------------------------------------------------------------------------------------------------------------------------





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
while getopts ":hs:b:f:d" OPT
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





#-------------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS
#-------------------------------------------------------------------------------------------------------------------------------

## Check if Sample Name variable exists
if [[ -z ${SAMPLE+x} ]] ## NOTE: ${VAR+x} is used for variable expansions, preventing unset variable error from set -o nounset. When $VAR is not set, we set it to "x" and throw the error.
then
        echo -e "$0 stopped at line ${LINENO}. \nREASON=Missing sample name option: -s"
        exit 1
fi

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.deliver_block_1.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.deliver_block_1.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
if [[ -z ${BAM+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing input BAM option: -b"
fi
if [[ ! -s ${BAM} ]]
then 
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Input sorted BAM file ${BAM} is empty or does not exist."
fi
if [[ ! -s ${BAM}.bai ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Sorted BAM index file ${BAM}.bai is empty or does not exist."
fi
if [[ -z ${JSON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing JSON option: -j"
fi
if [[ ! -s ${JSON} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Input JSON ${JSON} is empty or does not exist."
fi
if [[ -z ${DELIVERY_FOLDER+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing delivery folder option: -f"
fi
if [[ -d ${DELIVERY_FOLDER} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivery folder ${DELIVERY_FOLDER} already exists."
elif [[ -f ${DELIVERY_FOLDER} ]]
then 
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivery folder ${DELIVERY_FOLDER} is in fact a file."
fi

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## MAKE DELIVERY FOLDER
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Creating the Delivery folder."

## Copy the files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Creating Design Block 1 delivery folder. " ' INT TERM EXIT
mkdir -p ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Created the Design Block 1 delivery folder."







#-------------------------------------------------------------------------------------------------------------------------------
## DELIVERY
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Copying Design Block 1 outputs into Delivery folder."

## Copy the files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying BAM into delivery folder. " ' INT TERM EXIT
cp ${BAM} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Aligned sorted dedupped BAM delivered."


TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying BAM.BAI into delivery folder. " ' INT TERM EXIT
cp ${BAM}.bai ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Aligned sorted dedupped BAM.BAI delivered."

## Copy the JSON over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying JSON into delivery folder. " ' INT TERM EXIT
cp ${JSON} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Workflow JSON delivered."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index, and JSON. Open read permissions to the user group
BAM_NAME=`basename ${BAM}`
if [[ ! -s ${DELIVERY_FOLDER}/${BAM_NAME} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivered deduplicated BAM file ${DELIVERY_FOLDER}/${BAM_NAME} is empty."
fi
if [[ ! -s ${DELIVERY_FOLDER}/${BAM_NAME}.bai ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Deliveredt deduplicated BAM index file ${DELIVERY_FOLDER}/${BAM_NAME}.bai is empty."
fi
JSON_FILENAME=`basename ${JSON}` 
if [[ ! -s ${DELIVERY_FOLDER}/${JSON_FILENAME} ]]
then
       EXITCODE=1
       logError "$0 stopped at line ${LINENO}. \nREASON=Delivered workflow JSON file ${DELIVERY_FOLDER}/${JSON_FILENAME} is empty."
fi

chmod g+r ${DELIVERY_FOLDER}/${BAM_NAME}
chmod g+r ${DELIVERY_FOLDER}/${BAM_NAME}.bai
chmod g+r ${DELIVERY_FOLDER}/${JSON_FILENAME}


logInfo "[DELIVERY] Design Block 1 delivered. Have a nice day."

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
