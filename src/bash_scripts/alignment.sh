#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## alignment.sh MANIFEST, USAGE DOCS, SET CHECKS
#-------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST

*****************************************************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*****************************************************************************

MANIFEST
echo -e "\n${MANIFEST}"








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
                   -S           </path/to/sentieon> 
                   -L		<sentieon_license>
                   -t           <threads> 
                   -P		paired-end reads (true/false)
                   -d           turn on debug mode

 EXAMPLES:
 alignment.sh -h
 alignment.sh -g readgroup_ID -s sample -p platform -l read1.fq -r read2.fq -G reference.fa -S /path/to/sentieon_directory -L sentieon_license_number -t 12 -P true -d

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
while getopts ":hg:s:p:l:r:G:S:L:t:P:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage
                        echo -e "\n${DOCS}\n"
                        exit 0
			;;
                g )  # Read group ID. String variable invoked with -g
                        GROUP=${OPTARG}
                        ;;
                s )  # Sample name. String variable invoked with -s
                        SAMPLE=${OPTARG}
                        ;;
                p )  # Sequencing platform. String variable invoked with -p
                        PLATFORM=${OPTARG}
                        ;;
                l )  # Full path to input read 1. String variable invoked with -l
                        INPUT1=${OPTARG}
                        ;;
                r )  # Full path to input read 2. String variable invoked with -r
                        INPUT2=${OPTARG}
                        ;;
                G )  # Full path to referance genome fasta file. String variable invoked with -G
                        REFGEN=${OPTARG}
                        ;;
                S )  # Full path to sentieon directory. Invoked with -S
                        SENTIEON=${OPTARG}
                        ;;
		L )  # Sentieon license. Invoked with -L
			LICENSE=${OPTARG}
			;;
                t )  # Number of threads available. Integer invoked with -t
                        THR=${OPTARG}
                        ;;
                P )  # Is this a paired-end process? Boolean variable [true/false] invoked with -P
                        IS_PAIRED_END=${OPTARG}
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
ERRLOG=${SAMPLE}.alignment.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.alignment.log

if [[ ! -z ${ERRLOG+x} ]]
then
        echo -e "\nLog file ${ERRLOG} does not exist.\n"
        exit 1
fi

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"


## Check if input files, directories, and variables are non-zero
if [[ -z ${INPUT1+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing read 1 option: -l"
fi
if [[ ! -s ${INPUT1} ]]
then 
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input read 1 file ${INPUT1} is empty or does not exist."
fi
if [[ -z ${INPUT2+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing read 2 option: -r. If running a single-end job, set -r null in command."
fi
if [[ -z ${IS_PAIRED_END+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing paired-end option: -P"
fi
if [[ "${IS_PAIRED_END}" != true ]] && [[ "${IS_PAIRED_END}" != false ]]
then
	logError "$0 stopped at line ${LINENO}. \nREASON=Incorrect argument for paired-end option -P. Must be set to true or false."
fi
if [[ "${IS_PAIRED_END}" == true ]]
then
	if [[ ! -s ${INPUT2} ]]
	then
		EXITCODE=1
        	logError "$0 stopped at line $LINENO. \nREASON=Input read 2 file ${INPUT2} is empty or does not exist."
	fi
fi
if [[ -z ${REFGEN+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing reference genome option: -G"
fi
if [[ ! -s ${REFGEN} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REFGEN} is empty or does not exist."
fi
if [[ -z ${SENTIEON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Sentieon path option: -S"
fi
if [[ ! -d ${SENTIEON} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=BWA directory ${SENTIEON} is not a directory or does not exist."
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

## Set output file names
OUT=${SAMPLE}.sam
SORTBAM=${SAMPLE}.aligned.sorted.bam
SORTBAMIDX=${SAMPLE}.aligned.sorted.bam.bai

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## READ ALIGNMENT
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[BWA-MEM] START."

## BWA-MEM command, run for each read against a reference genome.
## Allocates all available threads to the process.
######## ASK ABOUT INTERLEAVED OPTION. NOTE: CAN ADD LANE TO RG OR REMOVE STRING
if [[ "${IS_PAIRED_END}" == false ]] # Align single read to reference genome
then
	export SENTIEON_LICENSE=${LICENSE}
	trap 'logError " $0 stopped at line ${LINENO}. Sentieon BWA-MEM error in read alignment. " ' INT TERM EXIT
	${SENTIEON}/bin/bwa mem -M -R "@RG\tID:$GROUP\tSM:${SAMPLE}\tPL:${PLATFORM}" -K 10000000 -t ${THR} ${REFGEN} ${INPUT1} > ${OUT}
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

	if [[ ${EXITCODE} -ne 0 ]]
        then
                logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        fi
else # Paired-end reads aligned
	export SENTIEON_LICENSE=${LICENSE}
	trap 'logError " $0 stopped at line ${LINENO}. Sentieon BWA-MEM error in read alignment. " ' INT TERM EXIT
	${SENTIEON}/bin/bwa mem -M -R "@RG\tID:$GROUP\tSM:${SAMPLE}\tPL:${PLATFORM}" -K 10000000 -t ${THR} ${REFGEN} ${INPUT1} ${INPUT2} > ${OUT}
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

        if [[ ${EXITCODE} -ne 0 ]]
        then
                logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
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
trap 'logError " $0 stopped at line ${LINENO}. Sentieon BAM conversion and sorting error. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon util sort -t ${THR} --sam2bam -i ${OUT} -o ${SORTBAM} >> ${SAMPLE}.alignment.log 2>&1
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi
logInfo "[SENTIEON] Converted output to BAM format and sorted."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check if BAM and index were created. Open read permissions to the user group
if [[ ! -s ${SORTBAM} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Output sorted BAM ${SORTBAM} is empty."
fi
if [[ ! -s ${SORTBAMIDX} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Output sorted BAM index ${SORTBAMIDX} is empty."
fi

chmod g+r ${OUT}
chmod g+r ${SORTBAM}
chmod g+r ${SORTBAMIDX}

logInfo "[BWA-MEM] Finished alignment. Aligned reads found in BAM format at ${SORTBAM}."

rm ${OUT}

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
