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
                   -K		<chunk_size_in_bases> 
                   -S           </path/to/sentieon> 
                   -t           <threads> 
                   -P		paired-end reads (true/false)
                   -e           </path/to/env_profile_file>
                   -d           turn on debug mode

 EXAMPLES:
 alignment.sh -h
 alignment.sh -g readgroup_ID -s sample -p platform -l read1.fq -r read2.fq -G reference.fa -K 10000000 -S /path/to/sentieon_directory -t 12 -P true -e /path/to/env_profile_file -d

 NOTE: To prevent different results due to thread count, set -K to 10000000 as recommended by the Sentieon manual.

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
while getopts ":hg:s:p:l:r:G:K:S:t:P:e:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage
                        echo -e "\n${DOCS}\n"
                        exit 0
			;;
                g )  # Read group ID
                        GROUP=${OPTARG}
			checkArg
                        ;;
                s )  # Sample name
                        SAMPLE=${OPTARG}
			checkArg
                        ;;
                p )  # Sequencing platform
                        PLATFORM=${OPTARG}
			checkArg
                        ;;
                l )  # Full path to input read 1
                        INPUT1=${OPTARG}
			checkArg
                        ;;
                r )  # Full path to input read 2
                        INPUT2=${OPTARG}
			checkArg
                        ;;
                G )  # Full path to referance genome fasta file
                        REFGEN=${OPTARG}
			checkArg
                        ;;
		K )  # Chunk size in bases (10000000 to prevent different results based on thread count)
			CHUNK_SIZE=${OPTARG}
			checkArg
			;;
                S )  # Full path to sentieon directory
                        SENTIEON=${OPTARG}
			checkArg
                        ;;
                t )  # Number of threads available
                        THR=${OPTARG}
			checkArg
                        ;;
                P )  # Is this a paired-end process? [true/false] Invoked with -P
                        IS_PAIRED_END=${OPTARG}
			checkArg
                        ;;
                e )  # Path to file with environmental profile variables
                        ENV_PROFILE=${OPTARG}
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
ERRLOG=${SAMPLE}.alignment.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.align_sentieon.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"


## source the file with environmental profile variables
if [[ ! -z ${ENV_PROFILE+x} ]]
then
        source ${ENV_PROFILE}
else
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing environmental profile option: -e"
fi

## Check if input files, directories, and variables are non-zero
if [[ -z ${INPUT1+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing read 1 option: -l"
fi
if [[ ! -s ${INPUT1} ]]
then 
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Input read 1 file ${INPUT1} is empty or does not exist."
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
	EXITCODE=1
	logError "$0 stopped at line ${LINENO}. \nREASON=Incorrect argument for paired-end option -P. Must be set to true or false."
fi
if [[ "${IS_PAIRED_END}" == true ]]
then
	if [[ ! -s ${INPUT2} ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}. \nREASON=Input read 2 file ${INPUT2} is empty or does not exist."
	fi
	if [[ "${INPUT2}" == null ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}/ \nREASON=User specified Paired End option -P, but set read 2 option -r to null."
	fi
fi
if [[ "${IS_PAIRED_END}" == false ]]
then
	if [[  "${INPUT2}" != null ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}/ \nREASON=User specified Single End option, but did not set read 2 option -r to null."
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
        logError "$0 stopped at line ${LINENO}. \nREASON=Reference genome file ${REFGEN} is empty or does not exist."
fi
if [[ -z ${CHUNK_SIZE+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line ${LINENO}. \nREASON=Missing read group option: -K\nSet -K 10000000 to prevent different results based on thread count."
fi
if [[ ${CHUNK_SIZE} != 10000000 ]]
then
	logWarn "[BWA-MEM] Chunk size option -K set to ${CHUNK_SIZE}. When this option is not set to 10000000, there may be different results per run based on different thread counts."
fi
if [[ -z ${GROUP+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing read group option: -g"
fi
if [[ -z ${PLATFORM+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing sequencing platform option: -p"
fi
if [[ -z ${SENTIEON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Sentieon path option: -S"
fi
if [[ ! -d ${SENTIEON} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=BWA directory ${SENTIEON} is not a directory or does not exist."
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
TOOL_LOG=${SAMPLE}.align_sentieon.log


#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## READ ALIGNMENT
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[BWA-MEM] START."

## BWA-MEM command, run for each read against a reference genome. 
if [[ "${IS_PAIRED_END}" == false ]] # Align single read to reference genome
then
	TRAP_LINE=$(($LINENO + 1))
	trap 'logError " $0 stopped at line ${TRAP_LINE}. Sentieon BWA-MEM error in read alignment. " ' INT TERM EXIT
	${SENTIEON}/bin/bwa mem -M -R "@RG\tID:$GROUP\tSM:${SAMPLE}\tPL:${PLATFORM}" -K ${CHUNK_SIZE} -t ${THR} ${REFGEN} ${INPUT1} > ${OUT} 2>>${TOOL_LOG}
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

	if [[ ${EXITCODE} -ne 0 ]]
        then
                logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        fi
else # Paired-end reads aligned
	TRAP_LINE=$(($LINENO + 1))
	trap 'logError " $0 stopped at line ${TRAP_LINE}. Sentieon BWA-MEM error in read alignment. " ' INT TERM EXIT
	${SENTIEON}/bin/bwa mem -M -R "@RG\tID:$GROUP\tSM:${SAMPLE}\tPL:${PLATFORM}" -K ${CHUNK_SIZE} -t ${THR} ${REFGEN} ${INPUT1} ${INPUT2} > ${OUT} 2>>${TOOL_LOG} 
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

        if [[ ${EXITCODE} -ne 0 ]]
        then
                logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
	fi
fi

if [[ ! -s ${OUT} ]]
then
	EXITCODE=1
	logError "$0 stopped at line ${LINENO}. \nREASON=Output SAM ${OUT} is empty."
fi


logInfo "[BWA-MEM] Aligned reads ${SAMPLE} to reference ${REFGEN}."

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## BAM CONVERSION AND SORTING
#-------------------------------------------------------------------------------------------------------------------------------

## Convert SAM to BAM and sort
logInfo "[SENTIEON] Converting SAM to BAM..."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Sentieon BAM conversion and sorting error. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon util sort -t ${THR} --sam2bam -i ${OUT} -o ${SORTBAM} >> ${TOOL_LOG} 2>&1
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
        logError "$0 stopped at line ${LINENO}. \nREASON=Output sorted BAM ${SORTBAM} is empty."
fi
if [[ ! -s ${SORTBAMIDX} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Output sorted BAM index ${SORTBAMIDX} is empty."
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
