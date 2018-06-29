#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## trim_sequences.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Trim input sequences using Cutadapt. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 trim_sequences.sh -s 		<sample_name> 
                   -A 		<adapters.fa> 
                   -l 		<read1.fq> 
                   -r 		<read2.fq> 
                   -C 		</path/to/cutadapt> 
                   -t 		<threads> 
                   -P 		paired-end reads (true/false)
                   -d 		turn on debug mode 

 EXAMPLES:
 trim_sequences.sh -h
 trim_sequences.sh -s sample -l read1.fq -r read2.fq -A adapters.fa -C /path/to/cutadapt_directory -t 12 -P true -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=trim_sequences.sh
SGE_JOB_ID=TBD   # placeholder until we parse job ID
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

while getopts ":hl:r:A:C:t:P:s:d" OPT
do
	case ${OPT} in
		h )  # Flag to display usage
			echo -e "\n${DOCS}\n"
			exit 0
			;;
		l )  # Full path to input read 1. String variable invoked with -l
			INPUT1=${OPTARG}
			;;
		r )  # Full path to input read 2. String variable invoked with -r
			INPUT2=${OPTARG}
			;;
		A )  # Full path to adapters fasta file. String variable invoked with -A
			ADAPTERS=${OPTARG}
			;;
		C )  # Full path to cutadapt installation directory. String variable invoked with -C
			CUTADAPT=${OPTARG}
			;;
		t )  # Number of threads available. Integer invoked with -t
			THR=${OPTARG}
			;;
		P )  # Is this a paired-end process? Boolean variable [true/false] invoked with -P
			IS_PAIRED_END=${OPTARG}
			;;
		s )  # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			;;
		d )  # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text
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

## Create log for JOB_ID/script and tool
ERRLOG=${SAMPLE}.trimming.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.cutadapt.log

## Send manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

if [[ -z ${ADAPTERS+x} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing adapters file option: -A"
fi
if [[ ! -s ${ADAPTERS} ]]  
then
	EXITCODE=1
	logError "$0 stopped at line ${LINENO}. \nREASON=Adapters fasta file ${ADAPTERS} is empty or does not exist."
fi
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
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing paired-end option: -P. If running a single-end job, set -P false in command."
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
fi
if [[ -z ${CUTADAPT+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing CutAdapt software path option: -C"
fi
if [[ ! -d ${CUTADAPT} ]]
then
	EXITCODE=1
	logError "$0 stopped at line ${LINENO}. \nREASON=Cutadapt directory ${CUTADAPT} is not a directory or does not exist."
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

## Parse filename without full path
OUT1=${SAMPLE}.read1.trimmed.fq.gz
if  [[ "${IS_PAIRED_END}" == false ]]  # If single-end, we do not need a second output trimmed read
then
	OUT2=null
else
	OUT2=${SAMPLE}.read2.trimmed.fq.gz
fi
#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## CUTADAPT READ TRIMMING
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[CUTADAPT] START."

## Cutadapt command, run for each fastq and each adapter sequence in the adapter FASTA file.
## Allocates half of the available threads to each process.
if [[ "${IS_PAIRED_END}" == false ]]  # if single-end reads file
then
	# Trim single-end reads

	trap 'logError " $0 stopped at line ${LINENO}. Cutadapt Read 1 failure. " ' INT TERM EXIT
	${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=${THR} -o ${OUT1} ${INPUT1} >> ${SAMPLE}.cutadapt.log 2>&1
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

	if [[ ${EXITCODE} -ne 0 ]]
	then
		logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
	fi
	logInfo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT log: ${SAMPLE}.cutadapt.log"
else
	# Trimming reads with Cutadapt in paired-end mode. -a and -A specify forward and reverse adapters, respectively. -p specifies output for read2 
	
	trap 'logError " $0 stopped at line ${LINENO}. Cutadapt Read 1 and 2 failure. " ' INT TERM EXIT
	${CUTADAPT}/cutadapt -a file:${ADAPTERS} -A file:${ADAPTERS} --cores=${THR} -p ${OUT2} -o ${OUT1} ${INPUT1} ${INPUT2} >> ${SAMPLE}.cutadapt.log 2>&1
	EXITCODE=$?
	trap - INT TERM EXIT

	if [[ ${EXITCODE} -ne 0 ]]
	then
		logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Cutadapt Read 1 and 2 failure."
	fi
	logInfo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT log: ${SAMPLE}.cutadapt.log"
fi


#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for file creation
if [[ ! -s ${OUT1} ]]
then
	EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Output trimmed read 1 file ${OUT1} is empty."
fi
if [[ "${IS_PAIRED_END}" == true ]]
then
	if [[ ! -s ${OUT2} ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}. \nREASON=Output trimmed read 2 file ${OUT2} is empty."
	fi
fi

## Open read permissions to the user group
if [[ "${IS_PAIRED_END}" == false ]]
then
	chmod g+r ${OUT1}
else
	chmod g+r ${OUT1}
	chmod g+r ${OUT2}
fi

logInfo "[CUTADAPT] Finished trimming adapter sequences."

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
