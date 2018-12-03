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
                   -e		</path/to/env_profile_file>
                   -F           </path/to/shared_functions.sh>
                   -d 		turn on debug mode 

 EXAMPLES:
 trim_sequences.sh -h
 trim_sequences.sh -s sample -l read1.fq -r read2.fq -A adapters.fa -C /path/to/cutadapt_directory -t 12 -P true -e /path/to/env_profile_file -F /path/to/shared_functions.sh -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=trim_sequences.sh
SGE_JOB_ID=TBD   # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID




#-------------------------------------------------------------------------------------------------------------------------------
## LOGGING FUNCTIONS
#-------------------------------------------------------------------------------------------------------------------------------

function checkArg()
{
    if [[ "${OPTARG}" == -* ]]; then
        echo -e "\nError with option -${OPT} in command. Option passed incorrectly or without argument.\n"
        echo -e "\n${DOCS}\n"
        exit 1;
    fi
}






#-------------------------------------------------------------------------------------------------------------------------------
## GETOPTS ARGUMENT PARSER
#-------------------------------------------------------------------------------------------------------------------------------

## Check if no arguments were passed
if (($# == 0))
then
	echo -e "\nNo arguments passed.\n\n${DOCS}\n"
	exit 1
fi

while getopts ":hl:r:A:C:t:P:s:e:F:d" OPT
do
	case ${OPT} in
		h )  # Flag to display usage
			echo -e "\n${DOCS}\n"
			exit 0
			;;
		l )  # Full path to input read 1
			INPUT1=${OPTARG}
			checkArg
			;;
		r )  # Full path to input read 2
			INPUT2=${OPTARG}
			checkArg
			;;
		A )  # Full path to adapters fasta file
			ADAPTERS=${OPTARG}
			checkArg
			;;
		C )  # Full path to cutadapt installation directory
			CUTADAPT=${OPTARG}
			checkArg
			;;
		t )  # Number of threads available
			THR=${OPTARG}
			checkArg
			;;
		P )  # Is this a paired-end process? [true/false] invoked with -P.
			IS_PAIRED_END=${OPTARG}
			checkArg
			;;
		s )  # Sample name
			SAMPLE=${OPTARG}
			checkArg
			;;
		e )  # Path to file with environmental profile variables
			ENV_PROFILE=${OPTARG}
			checkArg
			;;
                F )  # Path to shared_functions.sh
                        SHARED_FUNCTIONS=${OPTARG}
                        checkArg
                        ;;
		d )  # Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d.
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
## PRECHECK FOR INPUTS AND OPTIONS
#-------------------------------------------------------------------------------------------------------------------------------


source ${SHARED_FUNCTIONS}

## Check if Sample Name variable exists
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Create log for JOB_ID/script and tool
ERRLOG=${SAMPLE}.trimming.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.cutadapt.log

## Send manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with environmental profile variables
checkVar "${ENV_PROFILE+x}" "Missing environmental profile option: -e" $LINENO
source ${ENV_PROFILE}

##  Check if input files, directories, and variables are non-zero
checkVar "${ADAPTERS+x}" "Missing adapters file option: -A" $LINENO
checkFile ${ADAPTERS} "Adapters fasta file ${ADAPTERS} is empty or does not exist." $LINENO
checkVar "${INPUT1+x}" "Missing read 1 option: -l" $LINENO
checkFile ${INPUT1} "Input read 1 file ${INPUT1} is empty or does not exist." $LINENO
checkVar "${INPUT2+x}" "Missing read 2 option: -r. If running a single-end job, set -r null in command." $LINENO

checkVar "${IS_PAIRED_END+x}" "Missing paired-end option: -P" $LINENO

if [[ "${IS_PAIRED_END}" != true ]] && [[ "${IS_PAIRED_END}" != false ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Incorrect argument for paired-end option -P. Must be set to true or false."
fi
if [[ "${IS_PAIRED_END}" == true ]]
then
        checkFile ${INPUT2} "Input read 2 file ${INPUT2} is empty or does not exist. If running a single-end job, set -r null in command." $LINENO
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
		logError "$0 stopped at line ${LINENO}/ \nREASON=User specified Single End option, but did not set read 2 option -r to null."
	fi
fi

checkVar "${CUTADAPT+x}" "Missing CutAdapt software path option: -C" $LINENO
checkDir ${CUTADAPT} "Cutadapt directory ${CUTADAPT} is not a directory or does not exist." $LINENO
checkVar "${THR+x}" "Missing threads option: -t" $LINENO





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
	TRAP_LINE=$(($LINENO + 1))
	trap 'logError " $0 stopped at line ${TRAP_LINE}. Cutadapt Read 1 failure. " ' INT TERM EXIT
	${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=${THR} -o ${OUT1} ${INPUT1} >> ${SAMPLE}.cutadapt.log 2>&1
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

	checkExitcode ${EXITCODE} $LINENO
	logInfo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT log: ${SAMPLE}.cutadapt.log"
else
	# Trimming reads with Cutadapt in paired-end mode. -a and -A specify forward and reverse adapters, respectively. -p specifies output for read2 
	TRAP_LINE=$(($LINENO + 1))
	trap 'logError " $0 stopped at line ${TRAP_LINE}. Cutadapt Read 1 and 2 failure. " ' INT TERM EXIT
	${CUTADAPT}/cutadapt -a file:${ADAPTERS} -A file:${ADAPTERS} --cores=${THR} -p ${OUT2} -o ${OUT1} ${INPUT1} ${INPUT2} >> ${SAMPLE}.cutadapt.log 2>&1
	EXITCODE=$?
	trap - INT TERM EXIT

        checkExitcode ${EXITCODE} $LINENO
	logInfo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT log: ${SAMPLE}.cutadapt.log"
fi







#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for file creation
checkFile ${OUT1} "Output trimmed read 1 file ${OUT1} is empty." ${LINENO}
if [[ "${IS_PAIRED_END}" == true ]]
then
        checkFile ${OUT2} "Output trimmed read 2 file ${OUT2} is empty." ${LINENO}
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
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
