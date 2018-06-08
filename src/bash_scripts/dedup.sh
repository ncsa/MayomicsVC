#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## dedup.sh MANIFEST, USAGE DOCS, SET CHECKS
#-------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST
*******************************************
`readlink -m $0` was called by: `whoami` on `date`
command line input: ${@}
*******************************************
MANIFEST
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
                   -O           <output_directory> 
                   -S           </path/to/sentieon> 
                   -L		<sentieon_license>
                   -t           <threads> 
                   -e           </path/to/error_log> 
                   -d           debug_mode (true/false)

 EXAMPLES:
 dedup.sh -h
 dedup.sh -s sample -b aligned.sorted.bam -O /path/to/output_directory -S /path/to/sentieon_directory -L sentieon_license_number -t 12 -SE false -e /path/to/error.log -d true

#############################################################################

DOCS

set -o errexit
set -o pipefail
#set -o nounset

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

## Input and Output parameters
while getopts ":hs:b:O:S:L:t:e:d:" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
			echo " "
                        echo "Usage:"
			echo " "
                        echo "  bash dedup.sh -h       Display this help message."
                        echo "  bash dedup.sh [-s sample_name] [-b <aligned.sorted.bam>] [-O <output_directory>] [-S </path/to/sentieon>] [-L <sentieon_license>] [-t threads] [-e </path/to/error_log>] [-d debug_mode [false]]"
			echo " "
			exit 0;
                        ;;
		s )  # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			echo -e ${SAMPLE}
			;;
                b )  # Full path to the input BAM file. String variable invoked with -b
                        INPUTBAM=${OPTARG}
                        echo -e ${INPUTBAM}
                        ;;
                O )  # Output directory. String variable invoked with -O
                        OUTDIR=${OPTARG}
                        echo -e ${OUTDIR}
                        ;;
                S )  # Full path to sentieon directory. String variable invoked with -S
                        SENTIEON=${OPTARG}
                        echo -e ${SENTIEON}
                        ;;
		L )  # Sentieon license number. Invoked with -L
			LICENSE=${OPTARG}
			echo -e ${LICENSE}
			;;
                t )  # Number of threads available. Integer invoked with -t
                        THR=${OPTARG}
                        echo -e ${THR}
                        ;;
                e )  # Full path to error log file. String variable invoked with -e
                        ERRLOG=${OPTARG}
                        echo -e ${ERRLOG}
                        ;;
                d )  # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text
                        DEBUG=${OPTARG}
                        echo -e ${DEBUG}
                        ;;
        esac
done

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS
#-------------------------------------------------------------------------------------------------------------------------------

## Write manifest to log
echo "${MANIFEST}" >> ${ERRLOG}

## Turn on Debug Mode to print all code
if [[ ${DEBUG} == true ]]
then
	logInfo "Debug mode is ON."
        set -x
fi

## Check if input files, directories, and variables are non-zero
if [[ ! -d ${OUTDIR} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Output directory ${OUTDIR} does not exist."
        exit 1;
fi
if [[ ! -s ${INPUTBAM} ]]
then 
        logError "$0 stopped at line $LINENO. \nREASON=Input sorted BAM file ${INPUTBAM} is empty."
	exit 1;
fi
if [[ ! -s ${INPUTBAM}.bai ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Sorted BAM index file ${INPUTBAM}.bai is empty."
        exit 1;
fi
if [[ ! -d ${SENTIEON} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon directory ${SENTIEON} does not exist."
	exit 1;
fi
if (( ${THR} % 2 != 0 ))
then
	logWarn "Threads set to an odd integer. Subtracting 1 to allow for parallel, even threading."
	THR=$((THR-1))
fi

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Defining file names
samplename=${SAMPLE}
SCORETXT=${OUTDIR}/${SAMPLE}.score.txt
OUT=${OUTDIR}/${SAMPLE}.deduped.bam
OUTBAMIDX=${OUTDIR}/${SAMPLE}.deduped.bam.bai
DEDUPMETRICS=${OUTDIR}/${SAMPLE}.dedup_metrics.txt

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## DEDUPLICATION STAGE
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[SENTIEON] Collecting info to deduplicate BAM with Locus Collector."

## Locus Collector command
export SENTIEON_LICENSE=${LICENSE}
${SENTIEON}/bin/sentieon driver -t ${THR} -i ${INPUTBAM} --algo LocusCollector --fun score_info ${SCORETXT} 
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
	exit ${EXITCODE};
fi
logInfo "[SENTIEON] Locus Collector finished; starting Dedup."

## Dedup command (Note: optional --rmdup flag will remove duplicates; without, duplicates are marked but not removed)
export SENTIEON_LICENSE=${LICENSE}
${SENTIEON}/bin/sentieon driver -t ${THR} -i ${INPUTBAM} --algo Dedup --score_info ${SCORETXT} --metrics ${DEDUPMETRICS} ${OUT}
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi
logInfo "[SENTIEON] Deduplication Finished."
logInfo "[SENTIEON] Deduplicated BAM found at ${OUT}"

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index. Open read permissions to the user group
if [[ ! -s ${OUT} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Output deduplicated BAM file ${OUT} is empty."
        exit 1;
fi
if [[ ! -s ${OUTBAMIDX} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Output deduplicated BAM index file ${OUTBAMIDX} is empty."
        exit 1;
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
