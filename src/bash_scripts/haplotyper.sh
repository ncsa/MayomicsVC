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
		   -G	<reference_genome>
		   -t	<threads>
		   -b	<sorted.deduped.realigned.bam>
		   -D	<dbsnp.vcf>
		   -r	<recal_data.table>
		   -o	<extra_haplotyper_options>
                   -e   </path/to/env_profile_file>
		   -d   turn on debug mode

 EXAMPLES:
 Haplotyper.sh -h
 Haplotyper.sh -s sample -S /path/to/sentieon_directory -G reference.fa -t 12 -b sorted.deduped.realigned.recalibrated.bam -D dbsnp.vcf -r recal_data.table -o "'--emit_mode variant --gq_bands 1-60,60-99/19,99 --min_base_qual 10 --pcr_indel_model CONSERVATIVE --phasing 1 --ploidy 2 --prune_factor 2'" -e /path/to/env_profile_file -d 

NOTE: In order for getops to read in a string arguments for -o (extra_haplotyper_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.
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

    if [[ -z ${EXITCODE+x} ]]; then
        EXITCODE=1
    fi

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
        echo -e "\n${DOCS}\n "
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

while getopts ":hs:S:G:t:b:D:r:o:e:d" OPT
do
	case ${OPT} in 
		h ) # flag to display help message
			echo -e "\n${DOCS}\n "
			exit 0;
			;;
		s ) # Sample name
			SAMPLE=${OPTARG}
			checkArg
			;;
		S ) # Full path to sentieon directory
			SENTIEON=${OPTARG}
			checkArg
			;;
		G ) # Full path to referance genome fasta file
			REF=${OPTARG}
			checkArg
			;;
		t ) # Number of threads available
			NTHREADS=${OPTARG}
			checkArg
			;;
		b ) # Full path to BAM used as input
			INPUTBAM=${OPTARG}
			checkArg
			;;
		D ) # Full path to DBSNP file
			DBSNP=${OPTARG}
			checkArg
			;;
		r ) #Full path to the recal_data.table created in the BQSR step
			RECAL=${OPTARG}
			checkArg
			;;
		o ) #Extra options and arguments to haplotyper, input as a long string, can be empty if desired
			HAPLOTYPER_OPTIONS=${OPTARG}
			checkArg
			;;
                e )  # Path to file with environmental profile variables
                        ENV_PROFILE=${OPTARG}
                        checkArg
                        ;;
		d ) # Turn on debug mode. Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d
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

## source the file with environmental profile variables
if [[ ! -z ${ENV_PROFILE+x} ]]
then
        source ${ENV_PROFILE}
else
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing environmental profile option: -e"
fi


## Check if Sentieon path is present
if [[ -z ${SENTIEON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Sentieon path option: -S"
fi

## Check if the Sentieon executable is present.
if [[ ! -d ${SENTIEON} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon directory ${SENTIEON} is not a directory or does not exist."
fi

#Check if the number of threads is present.
if [[ -z ${NTHREADS+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing threads option: -t"
fi

## Check if the reference option was passed in.
if [[ -z ${REF+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing reference genome option: -G"
fi

## Check if the reference file exits
if [[ ! -s ${REF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REF} is empty or does not exist."
fi

## Check if the input BAM option was passed in.
if [[ -z ${INPUTBAM+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing input BAM option: -b"
fi

## Check if the BAM input file is present.
if [[ ! -s ${INPUTBAM} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input BAM ${INPUTBAM} is empty or does not exist."
fi

## Check if the input BAM index file is present
if [[ ! -s ${INPUTBAM}.bai ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input BAM index ${INPUTBAM} is empty or does not exist."
fi

## Check if the dbSNP option was passed in.
if [[ -z ${DBSNP+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing dbSNP option: -D"
fi


## Check if dbSNP file is present.
if [[ ! -s ${DBSNP} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=DBSNP ${DBSNP} is empty or does not exist."
fi

## Check if recal data table is option was passed in
if [[  -z ${RECAL+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing RECAL_DATA.TABLE option: -r"
fi

if [[ -z ${HAPLOTYPER_OPTIONS+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing extra haplotyper options option: -o"
fi

## Check if the Recal_data.table file produced in BQSR is present
if [[ ! -s ${RECAL} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=RECAL_DATA.TABLE ${RECAL} is empty or does not exist."
fi

#--------------------------------------------------------------------------------------------------------------------------------------------------
HAPLOTYPER_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${HAPLOTYPER_OPTIONS}`










#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform Haplotyper with Sentieon.
#--------------------------------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[Haplotyper] START."


#Execute Sentieon with the Haplotyper algorithm
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in Sentieon Haplotyper. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${RECAL} --algo Haplotyper ${HAPLOTYPER_OPTIONS_PARSED} -d ${DBSNP} ${SAMPLE}.vcf >> ${SAMPLE}.haplotype_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT

if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in Sentieon Haplotyper."
fi

logInfo "[Haplotyper] Finished running successfully. Output: ${SAMPLE}.vcf"
#------------------------------------------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for the creation of the output VCF file
if [[ ! -s ${SAMPLE}.vcf ]] 
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Output VCF is empty."
fi

## Open read permissions to the user group
chmod g+r ${SAMPLE}.vcf

## Check for the creation of the output VCF index file 
if [[ ! -s ${SAMPLE}.vcf.idx ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Output VCF index is empty."
fi


## Open read permissions to the user group
chmod g+r ${SAMPLE}.vcf.idx


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
