#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## dedup.sh MANIFEST, USAGE DOCS, SET CHECKS
#-------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST

*****************************************************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*****************************************************************************

MANIFEST
echo ""
echo "${MANIFEST}"







read -r -d '' DOCS << DOCS

#############################################################################
#
# Deduplicate BAMs with Sentieon. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 dedup.sh          -s           <sample_name> 
                   -b		<aligned.sorted.bam>
                   -S           </path/to/sentieon> 
                   -L		<sentieon_license>
                   -t           <threads> 
                   -d           turn on debug mode

 EXAMPLES:
 dedup.sh -h
 dedup.sh -s sample -b aligned.sorted.bam -S /path/to/sentieon_directory -L sentieon_license_number -t 12 -d

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
while getopts ":hs:b:S:L:t:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo "\n${DOCS}\n"
			exit 0
                        ;;
		s )  # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			;;
                b )  # Full path to the input BAM file. String variable invoked with -b
                        INPUTBAM=${OPTARG}
                        ;;
                S )  # Full path to sentieon directory. String variable invoked with -S
                        SENTIEON=${OPTARG}
                        ;;
		L )  # Sentieon license number. Invoked with -L
			LICENSE=${OPTARG}
			;;
                t )  # Number of threads available. Integer invoked with -t
                        THR=${OPTARG}
                        ;;
                d )  # Turn on debug mode. Initiates 'set -x' to print all text
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
if [[ -z ${SAMPLE+x} ]]
then
        echo -e "$0 stopped at line ${LINENO}. \nREASON=Missing sample name option: -s"
        exit 1
fi

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.${SGE_JOB_ID}.log

if [[ ! -f ${ERRLOG} ]]
then
        echo -e "\nLog file ${ERRLOG} is not a file.\n"
        exit 1
fi


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
        logError "$0 stopped at line $LINENO. \nREASON=Input sorted BAM file ${INPUTBAM} is empty or does not exist."
fi
if [[ ! -s ${INPUTBAM}.bai ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sorted BAM index file ${INPUTBAM}.bai is empty or does not exist."
fi
if [[ -z ${SENTIEON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Sentieon path option: -S"
fi
if [[ ! -d ${SENTIEON} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon directory ${SENTIEON} is empty or does not exist."
fi
if [[ -z ${LICENSE+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Sentieon license option: -L"
fi
if [[ -z ${THR+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing threads option: -t"
fi

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Defining file names
samplename=${SAMPLE}
SCORETXT=${SAMPLE}.deduped.score.txt
OUT=${SAMPLE}.aligned.sorted.deduped.bam
OUTBAMIDX=${SAMPLE}.aligned.sorted.deduped.bam.bai
DEDUPMETRICS=${SAMPLE}.dedup_metrics.txt

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## DEDUPLICATION STAGE
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[SENTIEON] Collecting info to deduplicate BAM with Locus Collector."

## Locus Collector command
export SENTIEON_LICENSE=${LICENSE}
${SENTIEON}/bin/sentieon driver -t ${THR} -i ${INPUTBAM} --algo LocusCollector --fun score_info ${SCORETXT} >> ${SAMPLE}.dedup.log 2>&1
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[SENTIEON] Locus Collector finished; starting Dedup."

## Dedup command (Note: optional --rmdup flag will remove duplicates; without, duplicates are marked but not removed)
export SENTIEON_LICENSE=${LICENSE}
${SENTIEON}/bin/sentieon driver -t ${THR} -i ${INPUTBAM} --algo Dedup --score_info ${SCORETXT} --metrics ${DEDUPMETRICS} ${OUT} >> ${SAMPLE}.dedup.log 2>&1
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[SENTIEON] Deduplication Finished. Deduplicated BAM found at ${OUT}"

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index. Open read permissions to the user group
if [[ ! -s ${OUT} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Output deduplicated BAM file ${OUT} is empty."
fi
if [[ ! -s ${OUTBAMIDX} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Output deduplicated BAM index file ${OUTBAMIDX} is empty."
fi

chmod g+r ${OUT}
chmod g+r ${OUTBAMIDX}
chmod g+r ${DEDUPMETRICS}

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
