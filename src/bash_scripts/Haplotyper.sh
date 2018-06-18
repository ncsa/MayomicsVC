#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------------
## Haplotyper.sh MANIFEST, USAGE DOCS, SET CHECKS
#------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST
*******************************************
`readlink -m $0` was called by: `whoami` on `date`
command line input: ${@}
*******************************************
MANIFEST
echo "${MANIFEST}"

read -r -d '' DOCS << DOCS


#######################################################################################################################################################
#
# Perform Sentieon's Haplotyper variant caller on the bam produced in the Deduplication stage of the Mayomics workflow
# BQSR must be run before to this stage to calculate the required modification of the quality scores.
# Step 2/3 in Single Sample Variant Calling.
#
########################################################################################################################################################

 USAGE:
 Haplotyper.sh     -s 	<sample_name>
 		   -O	</path/to/output_dir>
		   -S	</path/to/sentieon>
		   -L	<sentieon_license>
		   -G	</path/to/ref.fa>
		   -t	<threads>
		   -b	</path/to/deDuped.bam>
		   -D	</path/to/dbsnp.vcf>
		   -r	</path/to/recal_data.table>
		   -e	</path/to/error_log>
		   -d   debug_mode [false]

 EXAMPLES:
 Haplotyper.sh -h
 Haplotyper.sh -s sample -O output_dir -S sention -L sentieon_License -G ref.fa -t 12 -b sample.bam -D dbsnp.vcf -r recal_data.table -e error.log

##########################################################################################################################################################


DOCS

set -o errexit
set -o pipefail
#set -o nounset

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

    if [[ -n ${LOG_FILE-x} ]]; then
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
#--------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------
## GETOPS ARGUMENT PARSER
#--------------------------------------------------------------------------------------------------------------------------------
while getopts ":hs:O:S:L:G:t:b:D:r:e:d" OPT
do
	case ${OPT} in 
		h ) # flag to display help message
			echo " "
			echo "Usage:"
			echo " "
			echo "	bash Haplotyper.sh -h       Display this help message."
			echo "	bash Haplotyper.sh  [-s <sample_name>] [-O </path/to/output_dir>] [-S </path/to/sentieon>] [-G </path/to/ref.fa>] [-t <threads>] [-b </path/to/deDuped.bam>] [-D </path/to/dbsnp.vcf>] [-r </path/to/recal_data.table>] [-e </path/to/error_log>] [-d debug_mode [false]]"
			echo " "
			exit 0;
			;;
		s ) # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			;;
		O ) # Output directory. String variable invoked with -O
			OUTDIR=${OPTARG}
			;;
		S ) # Full path to Sentieon executable. String variable invoked with -S
			SENTIEON=${OPTARG}
			;;
		L ) # Sentieon license number. Invoked with -L 
			LICENSE=${OPTARG}
			;;
		G ) # Full path to reference fasta. String variable invoked with -r
			REF=${OPTARG}
			;;
		t ) # Number of threads available. Integer invoked with -t
			NTHREADS=${OPTARG}
			;;
		b ) # Full path to DeDuped BAM used as input. String variable invoked with -i
			INPUTBAM=${OPTARG}
			;;
		D ) # Full path to DBSNP file. String variable invoked with -D.
			DBSNP=${OPTARG}
			;;
		r ) #Full path to the recal_data.table created in the BQSR step
			RECAL=${OPTARG}
			;;
		e ) # Full path to error log file. String variable invoked with -e
			ERRLOG=${OPTARG}
			;;
		d ) # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text.
			DEBUG=${OPTARG}
			;;
	esac
done
#---------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS 
#---------------------------------------------------------------------------------------------------------------------------
## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"


## Turn on Debug Mode to print all code
if [[ ${DEBUG} == true ]]
then
        logInfo "Debug mode is ON."
        set -x
fi


## Check if error log file is present.
if [[ ! -f ${ERRLOG} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Error log file ${ERRLOG} does not exist."
        exit 1
fi

## Check if sample name is present.
if [[ -z ${SAMPLE} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=String for sample name is not present."
fi

## Check if the output directory exists.
if [[ ! -d ${OUTDIR} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Output directory ${OUTDIR} does not exist."
fi


## Check if the Sentieon executable is present.
if [[ ! -f ${SENTIEON} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon executable ${REF} is not present or does not exist."
fi

#Check if the number of threads is present.
if [[ -z ${NTHREADS} ]]
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
if [[ -z ${LICENSE} ]]
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
${SENTIEON} driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${RECAL} --algo Haplotyper -d ${DBSNP} ${OUTDIR}/${SAMPLE}.vcf
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in Sentieon Haplotyper."
        exit ${EXITCODE};
fi
#------------------------------------------------------------------------------------------------------------------------------------


## Open read permissions to the user group
chmod g+r -R ${OUTDIR}


logInfo "[Haplotyper] Finished running successfully. Output for ${SAMPLE} can be found in ${OUTDIR}/"
#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
