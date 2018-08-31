#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## deliver_block_2a.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Deliver results of Design Block 2a: vcf for snps and indels, their index files, and workflow JSON. 
# Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 deliver_block_2a.sh      -r           RecalibratedVcf 
                          -j           WorkflowJSONfile
                          -f           </path/to/delivery_folder>
                          -d           turn on debug mode

 EXAMPLES:
 deliver_block_2a.sh -h
 deliver_block_2a.sh -r Recalibrated.vcf -j Workflow.json -f /path/to/delivery_folder -d

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
while getopts ":hs:r:j:f:d" OPT
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
                r )  # Full path to the recalibrated VCF file
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
ERRLOG=${SAMPLE}.dedup.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.dedup_sentieon.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
if [[ -z ${VCF+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing VCF option: -r"
fi
if [[ ! -s ${VCF} ]]
then 
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Input VCF file ${VCF} is empty or does not exist."
fi
if [[ ! -s ${VCF}.idx ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Input VCF index file ${VCF}.idx is empty or does not exist."
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
trap 'logError " $0 stopped at line ${TRAP_LINE}. Creating Design Block 2a delivery folder. " ' INT TERM EXIT
mkdir -p ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Created the Design Block 2a delivery folder."








#-------------------------------------------------------------------------------------------------------------------------------
## DELIVERY
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[DELIVERY] Copying Design Block 2a outputs into Delivery folder."


## Copy the snp files over
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying VCF into delivery folder. " ' INT TERM EXIT
cp ${VCF} ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Recalibrated VCF delivered."


TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Copying VCF.IDX into delivery folder. " ' INT TERM EXIT
cp ${VCF}.idx ${DELIVERY_FOLDER}
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[DELIVERY] Recalibrated VCF.IDX delivered."


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

## Check for creation of output VCF and index, and JSON. Open read permissions to the user group
VCF_NAME=`basename ${VCF}`
if [[ ! -s ${DELIVERY_FOLDER}/${VCF_NAME} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivered recalibrated VCF file ${DELIVERY_FOLDER}/${VCF_NAME} is empty."
fi
if [[ ! -s ${DELIVERY_FOLDER}/${VCF_NAME}.idx ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Delivered recalibrated VCF index file ${DELIVERY_FOLDER}/${VCF_NAME}.idx is empty."
fi

JSON_FILENAME=`basename ${JSON}` 
if [[ ! -s ${DELIVERY_FOLDER}/${JSON_FILENAME} ]]
then
       EXITCODE=1
       logError "$0 stopped at line ${LINENO}. \nREASON=Delivered workflow JSON file ${DELIVERY_FOLDER}/${JSON_FILENAME} is empty."
fi


chmod g+r ${DELIVERY_FOLDER}/${VCF_NAME}
chmod g+r ${DELIVERY_FOLDER}/${VCF_NAME}.idx
chmod g+r ${DELIVERY_FOLDER}/${JSON_FILENAME}



logInfo "[DELIVERY] Design Block 2a delivered. Have a nice day."

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
