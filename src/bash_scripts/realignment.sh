#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## realignment.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Realign reads using Sentieon Realigner. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 realignment.sh    -s           <sample_name> 
                   -b           <deduped.bam>
                   -G		<reference_genome>
                   -k		<known_sites>
                   -O           <output_directory> 
                   -S           </path/to/sentieon> 
                   -L		<sentieon_license>
                   -t           <threads> 
                   -e           </path/to/error_log> 
                   -d           debug_mode (true/false)

 EXAMPLES:
 realignment.sh -h
 realignment.sh -s sample -b sorted.deduped.bam -G reference.fa -k indels.vcf,indels2.vcf,indels3.vcf -O /path/to/output_directory -S /path/to/sentieon_directory -L sentieon_license_number -t 12 -e /path/to/error.log -d true

#############################################################################

DOCS









set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=realignment.sh
SGE_JOB_ID=TBD  # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## LOGGING FUNCTIONS
#-------------------------------------------------------------------------------------------------------------------------------

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
while getopts ":hs:b:G:k:O:S:L:t:e:d:" OPT
do
        case ${OPT} in
                h )  # Flag to display usage
			echo " "
                        echo "Usage:"
			echo " "
                        echo "  bash realignment.sh -h       Display this help message."
                        echo "  bash realignment.sh [-s <sample_name>] [-b <input_bam>] [-G <reference_genome>] [-k <known_sites>] [-O <output_directory>] [-S </path/to/Sentieon>] [-L <sentieon_license>] [-t threads] [-e </path/to/error_log>] [-d debug_mode [false]]"
			echo " "
                        exit 0;
			;;
                s )  # Sample name. String variable invoked with -s
                        SAMPLE=${OPTARG}
                        echo -e ${SAMPLE}
                        ;;
		b )  # Full path to the input deduped BAM. String variable invoked with -b
			DEDUPEDBAM=${OPTARG}
			echo -e ${DEDUPEDBAM}
			;;
                G )  # Full path to referance genome fasta file. String variable invoked with -G
                        REFGEN=${OPTARG}
                        echo -e ${REFGEN}
                        ;;
		k )  # Full path to known sites file. String variable invoked with -k
			KNOWN=${OPTARG}
			echo -e ${KNOWN}
			;;
                O )  # Output directory. String variable invoked with -O
                        OUTDIR=${OPTARG}
                        echo -e ${OUTDIR}
                        ;;
                S )  # Full path to sentieon directory. Invoked with -S
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
echo "${MANIFEST}" >> "${ERRLOG}"

## Turn on Debug Mode to print all code
if [[ "${DEBUG}" == true ]]
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
if [[ ! -s ${DEDUPEDBAM} ]]
then
	logError "$0 stopped at line $LINENO. \nREASON=Deduped BAM ${DEDUPEDBAM} is empty."
	exit 1;
fi
if [[ ! -s ${DEDUPEDBAM}.bai ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Deduped BAM index ${DEDUPEDBAM} is empty."
        exit 1;
fi
if [[ ! -s ${REFGEN} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REFGEN} is empty."
        exit 1;
fi
#if [[ ! -s ${KNOWN} ]]
#then
#	logError "$0 stopped at line $LINENO. \nREASON=Known sites file ${KNOWN} is empty."
#	exit 1;
#fi
if [[ ! -d ${SENTIEON} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=BWA directory ${SENTIEON} does not exist."
	exit 1;
fi
#if (( ${THR} % 2 != 0 ))
#then
#	logWarn "Threads set to an odd integer. Subtracting 1 to allow for parallel, even threading."
#	THR=$((THR-1))
#fi

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME AND OPTION PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Parse known sites list of multiple files. Create multiple -k flags for sentieon
SPLITKNOWN=`sed -e 's/,/ -k /g' <<< ${KNOWN}`
echo ${SPLITKNOWN}

## Parse filenames without full path
OUT=${OUTDIR}/${SAMPLE}.realigned.bam

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## REALIGNMENT STAGE
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[Realigner] START. Realigning deduped BAM. Using known sites at ${KNOWN}."

## Sentieon Realigner command.
## Allocates all available threads to the process.
export SENTIEON_LICENSE=${LICENSE}
${SENTIEON}/bin/sentieon driver -t ${THR} -r ${REFGEN} -i ${DEDUPEDBAM} --algo Realigner -k ${SPLITKNOWN} ${OUT}
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi
logInfo "[Realigner] Realigned reads ${SAMPLE} to reference ${REFGEN}. Realigned BAM located at ${OUT}."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of realigned BAM and index. Open read permissions to the user group
if [[ ! -s ${OUT} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Realigned BAM ${OUT} is empty."
        exit 1;
fi
if [[ ! -s ${OUT}.bai ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Realigned BAM ${OUT}.bai is empty."
        exit 1;
fi
chmod g+r ${OUT}

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
