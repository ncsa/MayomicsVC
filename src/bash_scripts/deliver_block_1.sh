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
                          -f           </path/to/delivery_folder>
                          -d           turn on debug mode

 EXAMPLES:
 deliver_block_1.sh -h
 deliver_block_1.sh -b aligned.sorted.dedupped.bam -f /path/to/delivery_folder -d

#############################################################################

DOCS








set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=dedup.sh
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
while getopts ":hb:f:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo -e "\n${DOCS}\n"
			exit 0
                        ;;
                b )  # Full path to the input BAM file
                        INPUTBAM=${OPTARG}
			checkArg
                        ;;
                e )  # Path to delivery folder
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

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.dedup.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.dedup_sentieon.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
if [[ -z ${INPUTBAM+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing input BAM option: -b"
fi
if [[ ! -s ${INPUTBAM} ]]
then 
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Input sorted BAM file ${INPUTBAM} is empty or does not exist."
fi
if [[ ! -s ${INPUTBAM}.bai ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Sorted BAM index file ${INPUTBAM}.bai is empty or does not exist."
fi
if [[ -z ${DELIVERY_FOLDER+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing delivery folder option: -f"
fi
if [[ ! -d ${DELIVERY_FOLDER} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivery folder ${DELIVERY_FOLDER} is not a directory or does not exist."
fi

#-------------------------------------------------------------------------------------------------------------------------------








#-------------------------------------------------------------------------------------------------------------------------------
## DELIVERY
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[SENTIEON] Collecting info to deduplicate BAM with Locus Collector."

## Copy the files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying BAM into delivery folder. " ' INT TERM EXIT
cp ${INPUTBAM} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Aligned sorted dedupped BAM delivered."


TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying BAM.BAI into delivery folder. " ' INT TERM EXIT
cp ${INPUTBAM}.bai ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Aligned sorted dedupped BAM.BAI delivered."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index. Open read permissions to the user group
if [[ ! -s ${DELIVERY_FOLDER}/${INPUTBAM} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivered deduplicated BAM file ${DELIVERY_FOLDER}/${INPUTBAM} is empty."
fi
if [[ ! -s ${DELIVERY_FOLDER}/${INPUTBAM}.bai ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Deliveredt deduplicated BAM index file ${DELIVERY_FOLDER}/${INPUTBAM}.bai is empty."
fi

chmod g+r ${DELIVERY_FOLDER}/${INPUTBAM}
chmod g+r ${DELIVERY_FOLDER}/${INPUTBAM}.bai

logInfo "[DELIVERY] Design Block 1 delivered. Have a nice day."

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
