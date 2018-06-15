#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## alignment.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Align sequences using Sentieon/BWA-MEM. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 alignment.sh      -g		<readgroup_ID>
                   -s           <sample_name> 
                   -p		<platform>
                   -l           <read1.fq> 
                   -r           <read2.fq>
                   -G		<reference_genome> 
                   -O           <output_directory> 
                   -S           </path/to/sentieon> 
                   -L		<sentieon_license>
                   -t           <threads> 
                   -P		single-end read (true/false)
                   -e           </path/to/error_log> 
                   -d           debug_mode (true/false)

 EXAMPLES:
 alignment.sh -h
 alignment.sh -g readgroup_ID -s sample -p platform -l read1.fq -r read2.fq -G reference.fa -O /path/to/output_directory -S /path/to/sentieon_directory -L sentieon_license_number -t 12 -P false -e /path/to/error.log -d true

#############################################################################

DOCS







set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=alignment.sh
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
while getopts ":hg:s:p:l:r:G:O:S:L:t:P:e:d:" OPT
do
        case ${OPT} in
                h )  # Flag to display usage
			echo " "
                        echo "Usage:"
			echo " "
                        echo "  bash alignment.sh -h       Display this help message."
                        echo "  bash alignment.sh [-g <readgroup_ID>] [-s <sample_name>] [-p <platform>] [-l <read1.fq>] [-r <read2.fq>] [-G <reference_genome>] [-O <output_directory>] [-S </path/to/Sentieon>] [-L <sentieon_license>] [-t threads] [-P single-end? (true/false)] [-e </path/to/error_log>] [-d debug_mode [false]]"
			echo " "
                        exit 0;
			;;
                g )  # Read group ID. String variable invoked with -g
                        GROUP=${OPTARG}
                        echo -e ${GROUP}
                        ;;
                s )  # Sample name. String variable invoked with -s
                        SAMPLE=${OPTARG}
                        echo -e ${SAMPLE}
                        ;;
                p )  # Sequencing platform. String variable invoked with -p
                        PLATFORM=${OPTARG}
                        echo -e ${PLATFORM}
                        ;;
                l )  # Full path to input read 1. String variable invoked with -l
                        INPUT1=${OPTARG}
                        echo -e ${INPUT1}
                        ;;
                r )  # Full path to input read 2. String variable invoked with -r
                        INPUT2=${OPTARG}
                        echo -e ${INPUT2}
                        ;;
                G )  # Full path to referance genome fasta file. String variable invoked with -G
                        REFGEN=${OPTARG}
                        echo -e ${REFGEN}
                        ;;
                O )  # Output directory. String variable invoked with -O
                        OUTDIR=${OPTARG}
                        echo -e ${OUTDIR}
                        ;;
                S )  # Full path to sentieon directory. Invoked with -S
                        SENTIEON=${OPTARG}
                        echo -e ${SENTIEON}
                        ;;
		L )  # Sentieon license. Invoked with -L
			LICENSE=${OPTARG}
			echo -e ${SENTIEON}
			;;
                t )  # Number of threads available. Integer invoked with -t
                        THR=${OPTARG}
                        echo -e ${THR}
                        ;;
                P )  # Is this a single-end process? Boolean variable [true/false] invoked with -P
                        IS_SINGLE_END=${OPTARG}
                        echo -e ${IS_SINGLE_END}
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
if [[ ! -s ${INPUT1} ]]
then 
        logError "$0 stopped at line $LINENO. \nREASON=Input read 1 file ${INPUT1} is empty."
	exit 1;
fi
if [[ ${IS_SINGLE_END} == false ]]
then
	if [[ ! -s ${INPUT2} ]]
	then
        	logError "$0 stopped at line $LINENO. \nREASON=Input read 2 file ${INPUT2} is empty."
		exit 1;
	fi
fi
if [[ ! -s ${REFGEN} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REFGEN} is empty."
        exit 1;
fi
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
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Set output file names
OUT=${OUTDIR}/${SAMPLE}.sam
SORTBAM=${OUTDIR}/${SAMPLE}.sorted.bam
SORTBAMIDX=${OUTDIR}/${SAMPLE}.sorted.bam.bai

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## READ ALIGNMENT
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[BWA-MEM] START."

## BWA-MEM command, run for each read against a reference genome.
## Allocates all available threads to the process.
######## ASK ABOUT INTERLEAVED OPTION. NOTE: CAN ADD LANE TO RG OR REMOVE STRING
if [[ "${IS_SINGLE_END}" == true ]]
then
	export SENTIEON_LICENSE=${LICENSE}
	${SENTIEON}/bin/bwa mem -M -R "@RG\tID:$GROUP\tSM:${SAMPLE}\tPL:${PLATFORM}" -K 100000000 -t ${THR} ${REFGEN} ${INPUT1} > ${OUT}
	EXITCODE=$?  # Capture exit code
	if [[ ${EXITCODE} -ne 0 ]]
        then
                logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
                exit ${EXITCODE};
        fi
else
	export SENTIEON_LICENSE=${LICENSE}
	${SENTIEON}/bin/bwa mem -M -R "@RG\tID:$GROUP\tSM:${SAMPLE}\tPL:${PLATFORM}" -K 100000000 -t ${THR} ${REFGEN} ${INPUT1} ${INPUT2} > ${OUT}
	EXITCODE=$?  # Capture exit code
        if [[ ${EXITCODE} -ne 0 ]]
        then
                logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
                exit ${EXITCODE};
	fi
fi
logInfo "[BWA-MEM] Aligned reads ${SAMPLE} to reference ${REFGEN}."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## BAM CONVERSION AND SORTING
#-------------------------------------------------------------------------------------------------------------------------------

## Convert SAM to BAM and sort
logInfo "[SENTIEON] Converting SAM to BAM..."
export SENTIEON_LICENSE=${LICENSE}
${SENTIEON}/bin/sentieon util sort -t ${THR} --sam2bam -i ${OUT} -o ${SORTBAM} >> ${ERRLOG}
EXITCODE=$?  # Capture exit code
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
	exit ${EXITCODE};
fi
logInfo "[SENTIEON] Converted output to BAM format and sorted."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check if BAM and index were created. Open read permissions to the user group
if [[ ! -s ${SORTBAM} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Output sorted BAM ${SORTBAM} is empty."
        exit 1;
fi
if [[ ! -s ${SORTBAMIDX} ]]
then
        logError "$0 stopped at line $LINENO. \nREASON=Output sorted BAM index ${SORTBAMIDX} is empty."
        exit 1;
fi

chmod g+r ${OUT}
chmod g+r ${SORTBAM}
chmod g+r ${SORTBAMIDX}

logInfo "[BWA-MEM] Finished alignment. Aligned reads found in BAM format at ${SORTBAM}."

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
