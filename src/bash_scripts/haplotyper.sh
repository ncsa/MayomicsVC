#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------------
## haplotyper.sh MANIFEST, USAGE DOCS, SET CHECKS
#------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST
*******************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*******************************************
MANIFEST
echo -e "\n${MANIFEST}"






read -r -d '' DOCS << DOCS


#######################################################################################################################################################
#
# Perform Sentieon's Haplotyper variant caller on the bam produced in the Deduplication stage of the Mayomics workflow.
# bqsr.sh must be run before to this stage to calculate the required modification of the quality scores.
# Step 2/3 in Single Sample Variant Calling.
#
########################################################################################################################################################

 USAGE:
 Haplotyper.sh     -s 	<sample_name>
		   -S	</path/to/sentieon>
		   -L	<sentieon_license>
		   -G	<reference_genome>
		   -t	<threads>
		   -b	<sorted.deduped.realigned.bam>
		   -D	<dbsnp.vcf>
		   -r	<recal_data.table>
		   -d   turn on debug mode

 EXAMPLES:
 Haplotyper.sh -h
 Haplotyper.sh -s sample -S sention -L sentieon_license -G ref.fa -t 12 -b sorted.deduped.realigned.bam -D dbsnp.vcf -r recal_data.table -d 

##########################################################################################################################################################


DOCS

set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=Haplotyper.sh
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


#SEVERITY=ERROR
function logError()
{
    local LEVEL="ERROR"
    local CODE="-1"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    >&2 _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"
    exit ${EXITCODE};
}


#SEVERITY=WARN
function logWarn()
{
    local LEVEL="WARN"
    local CODE="0"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"
}

#SEVERITY=INFO
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
        exit 1;
    fi
}

#--------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------
## GETOPS ARGUMENT PARSER
#--------------------------------------------------------------------------------------------------------------------------------

## Check if no arguments were passed
if (($# == 0))
then
        echo -e "\nNo arguments passed.\n\n${DOCS}\n"
        exit 1
fi

while getopts ":hs:S:L:G:t:b:D:r:d" OPT
do
	case ${OPT} in 
		h ) # flag to display help message
			echo -e "\n${DOCS}\n "
			exit 0;
			;;
		s ) # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			checkArg
			;;
		S ) # Full path to Sentieon. String variable invoked with -S
			SENTIEON=${OPTARG}
			checkArg
			;;
		L ) # Sentieon license number. Invoked with -L 
			LICENSE=${OPTARG}
			checkArg
			;;
		G ) # Full path to reference fasta. String variable invoked with -r
			REF=${OPTARG}
			checkArg
			;;
		t ) # Number of threads available. Integer invoked with -t
			NTHREADS=${OPTARG}
			checkArg
			;;
		b ) # Full path to DeDuped BAM used as input. String variable invoked with -i
			INPUTBAM=${OPTARG}
			checkArg
			;;
		D ) # Full path to DBSNP file. String variable invoked with -D.
			DBSNP=${OPTARG}
			checkArg
			;;
		r ) #Full path to the recal_data.table created in the BQSR step
			RECAL=${OPTARG}
			checkArg
			;;
		d ) # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text.
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
#---------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS 
#---------------------------------------------------------------------------------------------------------------------------
## Check if sample name is set
if [[ -z ${SAMPLE+x} ]] ## NOTE: ${VAR+x} is used for variable expansions, preventing unset variable error from set -o nounset. When $VAR is not set, we set it to "x" and throw the error.
then
	echo -e "$0 stopped at line $LINENO. \nREASON=Missing sample name option: -s"
	exit 1
fi

## Send Manifest to log
ERRLOG=${SAMPLE}.haplotyper.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.haplotype_sentieon.log

echo "${MANIFEST}" >> "${ERRLOG}"

## Check if the Sentieon executable is present.
if [[ ! -d ${SENTIEON} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon executable ${REF} is not present or does not exist."
fi

#Check if the number of threads is present.
if [[ -z ${NTHREADS+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Number of threads is not specified."
fi

## Check if the reference fasta file is present.
if [[ ! -f ${REF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome fasta file ${REF} is not present or does not exist."
fi

## Check if the BAM input file is present.
if [[ ! -f ${INPUTBAM} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input BAM ${INPUTBAM} is not present or does not exist."
fi

## Check if dbSNP file is present.
if [[ ! -f ${DBSNP} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=DBSNP ${DBSNP} is not present or does not exist."
fi

## Check if the Recal_data.table file produced in BQSR is present
if [[ ! -f ${RECAL} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=RECAL_DATA.TABLE ${RECAL} is not present or does not exist."
fi

## Check if Sentieon license string is present.
if [[ -z ${LICENSE+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon license ${LICENSE} is not present or does not exist."
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform Haplotyper with Sentieon.
#--------------------------------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[Haplotyper] START."

export SENTIEON_LICENSE=${LICENSE}

#Execute sentieon with the Haplotyper algorithm
trap 'logError " $0 stopped at line ${LINENO}. Cutadapt Read 1 failure. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${RECAL} --algo Haplotyper -d ${DBSNP} ${SAMPLE}.vcf >> ${SAMPLE}.haplotype_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in Sentieon Haplotyper."
        exit ${EXITCODE};
fi
#------------------------------------------------------------------------------------------------------------------------------------


## Open read permissions to the user group
chmod g+r ${SAMPLE}.vcf


logInfo "[Haplotyper] Finished running successfully. Output: ${SAMPLE}.vcf"
#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
