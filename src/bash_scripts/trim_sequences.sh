#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## trim_sequences.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Trim input sequences using Cutadapt. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 trim_sequences.sh -s 		<sample_name> 
                   -A 		<adapters.fa> 
                   -r 		<read1.fq> 
                   -R 		<read2.fq> 
                   -O 		<output_directory> 
                   -C 		</path/to/cutadapt> 
                   -t 		<threads> 
                   -SE 		single-end read (true/false)
                   -e 		</path/to/error_log> 
                   -d 		debug_mode (true/false)

 EXAMPLES:
 trim_sequences.sh -h
 trim_sequences.sh -s sample -r read1.fq -R read2.fq -A adapters.fa -O /path/to/output_directory -C /path/to/cutadapt_directory -t 12 -SE false -e /path/to/error.log -d true

#############################################################################

DOCS

set -o errexit
set -o pipefail
#set -o nounset

SCRIPT_NAME=trim_sequences.sh
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

while getopts ":hs:r:R:A:O:C:t:SE:e:d:" OPT
do
	case ${OPT} in
		h )  # Flag to display usage
			echo " "
			echo "Usage:"
			echo " "
			echo "	bash trim_sequences.sh -h	Display this help message."
			echo "	bash trim_sequences.sh [-s sample_name] [-r <read1.fq>] [-R <read2.fq>] [-A <adapters.fa>] [-O </path/to/output_directory>] [-C </path/to/cutadapt_directory>] [-t threads] [-SE single_end? (true/false)] [-e <error_log>] [-d debug_mode [false]]"
			echo " "
			exit 0;
			;;
		s )  # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			echo ${SAMPLE}
			;;
		r )  # Full path to input read 1. String variable invoked with -r
			INPUT1=${OPTARG}
			echo ${INPUT1}
			;;
		R )  # Full path to input read 2. String variable invoked with -r
			INPUT2=${OPTARG}
			echo ${INPUT2}
			;;
		A )  # Full path to adapters fasta file. String variable invoked with -A
			ADAPTERS=${OPTARG}
			echo ${ADAPTERS}
			;;
		O )  # Output directory. String variable invoked with -O
			OUTDIR=${OPTARG}
			echo ${OUTDIR}
			;;
		C )  # Full path to cutadapt directory. String variable invoked with -C
			CUTADAPT=${OPTARG}
			echo ${CUTADAPT}
			;;
		t )  # Number of threads available. Integer invoked with -t
			THR=${OPTARG}
			echo ${THR}
			;;
		SE )  # Is this a single-end process? Boolean variable [true/false] invoked with -SE
			IS_SINGLE_END=${OPTARG}
			echo ${IS_SINGLE_END}
			;;
		e )  # Full path to error log file. String variable invoked with -e
			ERRLOG=${OPTARG}
			echo ${ERRLOG}
			;;
		d )  # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text
			DEBUG=${OPTARG}
			echo ${DEBUG}
			;;
	esac
done

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS
#-------------------------------------------------------------------------------------------------------------------------------

## Turn on Debug Mode to print all code
if [[ ${DEBUG} == true ]]
then
	logInfo "Debug mode is ON."
	set -x
fi

## Check if input files, directories, and variables are non-zero
if [[ ! -d ${OUTDIR} ]]
then
        logError "$0 stopped at line ${LINENO}. \nREASON=Output directory ${OUTDIR} does not exist."
        exit 1;
fi
if [[ ! -s ${ADAPTERS} ]]  
then
	logError "$0 stopped at line ${LINENO}. \nREASON=Adapters fasta file ${ADAPTERS} is empty."
	exit 1;
fi
if [[ ! -s ${INPUT1} ]]  
then
	logError "$0 stopped at line ${LINENO}. \nREASON=Input read 1 file ${INPUT1} is empty."
	exit 1;
fi
if [[ ${IS_SINGLE_END} == false ]]
then
        if [[ ! -s ${INPUT2} ]]
        then
                logError "$0 stopped at line ${LINENO}. \nREASON=Input read 2 file ${INPUT2} is empty."
                exit 1;
        fi
fi
if [[ ! -d ${CUTADAPT} ]]
then
	logError "$0 stopped at line ${LINENO}. \nREASON=Cutadapt directory ${CUTADAPT} does not exist."
	exit 1;
fi
if (( ${THR} % 2 != 0 ))  ## This is checking if the number of threads is an odd number. If that is the case, we subtract 1 from the integer so the parallel processes can run on equal threads.
then
	logWarn "Threads set to an odd integer. Subtracting 1 to allow for parallel, even threading."
	THR=$((THR-1))
fi

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Parse filename without full path
full1=$INPUT1
full2=$INPUT2
READ1=${full1##*/} # Remove path from variable
READ2=${full2##*/}
read1=${READ1%%.*} # Remove all instances of file extensions
read2=${READ2%%.*}
OUT=${OUTDIR}/${SAMPLE}.trimmed.fq.gz
OUT1=${OUTDIR}/${SAMPLE}.read1.trimmed.fq.gz
OUT2=${OUTDIR}/${SAMPLE}.read2.trimmed.fq.gz

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## CUTADAPT READ TRIMMING
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[CUTADAPT] START."

## Cutadapt command, run for each fastq and each adapter sequence in the adapter FASTA file.
## Allocates half of the available threads to each process.
if [[ ${IS_SINGLE_END} = true ]]  # if single-end reads file
then
	# Trim reads
	${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=${THR} -o ${OUT} ${INPUT1} >> ${SAMPLE}.cutadapt.log 
	EXITCODE=$?  # Capture exit code
	if [[ ${EXITCODE} -ne 0 ]]
	then
		logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
		exit ${EXITCODE};
	fi
else 
	if [[ $THR = 0 ]] # if threads equals zero
	then
		${CUTADAPT}/cutadapt -a file:${ADAPTERS} -o ${OUT1} ${INPUT1} >> ${read1}.cutadapt.log &
		${CUTADAPT}/cutadapt -a file:${ADAPTERS} -o ${OUT2} ${INPUT2} >> ${read2}.cutadapt.log &
		wait
		EXITCODE=$?
                if [[ ${EXITCODE} -ne 0 ]]
                then
                        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
                        exit ${EXITCODE};
                fi
		
	else  # If threads does not equal zero

		${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=$((THR/2)) -o ${OUT1} ${INPUT1} >> ${read1}.cutadapt.log &
		${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=$((THR/2)) -o ${OUT2} ${INPUT2} >> ${read2}.cutadapt.log &
		wait
		EXITCODE=$?
                if [[ ${EXITCODE} -ne 0 ]]
                then
                        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
                        exit ${EXITCODE};
                fi
		
	fi
fi

logInfo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT logs: ${OUT}/${read1}.cutadapt.log ${OUT}/${read2}.cutadapt.log"

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Open read permissions to the user group
if [[ ${IS_SINGLE_END} = true ]]
then
	chmod g+r ${OUT}
else
	chmod g+r ${OUT1}
	chmod g+r ${OUT2}
fi

logInfo "[CUTADAPT] Finished trimming adapter sequences. Trimmed reads found at ${OUTDIR}/"

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
