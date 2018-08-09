#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## vqsr.sh MANIFEST, USAGE DOCS, SET CHECKS
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

#################################################################################################################################
#
# Perform VQSR (Variant Quality Score Recalibration on the VCF produced from the Haplotyper stage of the Mayomics Workflow.
# Step 3/3 in Single Sample Variant Calling.
#
#################################################################################################################################

 USAGE:
 vqsr.sh -s <sample_name>
 	 -S </path/to/sentieon>
	 -L <Sentieon_license>
	 -G <reference_genome>
	 -V <sample.vcf>
	 -r <resource_string_for_SNPs>
	 -R <resource_string_for_INDELS>
	 -d turn on debug mode

 EXAMPLES:
 vqsr.sh -h
 vqsr.sh -s sample -S /path/to/sentieon_directory -L sentieon_license_number -G reference.fa -V sample.vcf -r "--resource 1000G.vcf --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 --resource omni.vcf --resource_param omni,known=false,training=true,truth=false,prior=12.0 --resource dbSNP.vcf --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 --resource hapmap.vcf --resource_param hapmap,known=false,training=true,truth=true,prior=15.0" -R "--resource mills.vcf --resource_param Mills,known=false,training=true,truth=true,prior=12.0 --resource dbSNP.vcf --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0" -d

##################################################################################################################################

DOCS











set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=vqsr.sh
SGE_JOB_ID=TBD  # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID

#------------------------------------------------------------------------------------------------------------------------------




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
while getopts ":hs:S:L:G:V:r:R:d" OPT
do
	case ${OPT} in
		h ) # Flag to display usage
			echo -e "\n${DOCS}\n"
			exit 0
			;;
		s ) # Sample name
			SAMPLE=${OPTARG}
			checkArg
			;;
		S ) # Full path to sentieon directory
			SENTIEON=${OPTARG}
			checkArg
			;;
		L ) # Sentieon license
			LICENSE=${OPTARG}
			checkArg
			;;
		G ) # Reference genome
			REF=${OPTARG}
			checkArg
			;;
		V ) # Sample VCF file output from Haplotyper
			SAMPLEVCF=${OPTARG}
			checkArg
			;;
		r ) # Resource string for SNPs
			RESOURCE_SNPS=${OPTARG}
			checkArg
			;;
		R ) # Resource string for INDELs
			RESOURCE_INDELS=${OPTARG}
			checkArg
			;;
		d ) # Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d.
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
## Check if sample name is present.
if [[ -z ${SAMPLE+x} ]] ## NOTE: ${VAR+x} is used for variable expansions, preventing unset variable error from set -o nounset. When $VAR is not set, we set it to "x" and throw the error.
then
	echo -e "$0 stopped at line ${LINENO}. \nREASON=Missing sample name option: -s"
	exit 1
fi

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.vqsr.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.vqsr_sentieon.log

## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if the Sentieon executable option was passed in
if [[ -z ${SENTIEON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon executable option: -S"
fi

## Check if the Sentieon executable is present
if [[ ! -d ${SENTIEON} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon directory ${SENTIEON} is not a directory or does not exist."
fi

## Check if the Sentieon license option was passed in
if [[ -z ${LICENSE+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon license option: -S"
fi

## Check if the reference option was passed in
if [[ -z ${REF+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing reference fasta file option: -G"
fi

## Check if the reference fasta file is present
if [[ ! -s ${REF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REF} is empty or does not exist."
fi

## Check if the sample VCF input file option was passed in
if [[ -z ${SAMPLEVCF+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing sample input VCF option: -V"
fi

## Check if the sample VCF input file is present
if [[ ! -s ${SAMPLEVCF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input sample vcf ${SAMPLEVCF} is empty or does not exist."
fi

## Check if the resource string for SNPs was passed in
if [[ -z ${RESOURCE_SNPS+x} ]] 
then 
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing resource string for SNPs option: -r"
fi

## Check if the resource string for INDELS was passed in
if [[ -z ${RESOURCE_INDELS+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing resource string for INDELs option: -R"
fi

#-------------------------------------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------------------
# VQSR STAGE
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[VQSR] START. Performing VQSR on VCF output from Haplotyper."

export SENTIEON_LICENSE=${LICENSE}


## Create the ANNOTATION argument
ANNOTATE_TEXT="--annotation DP --annotation QD --annotation FS --annotation SOR --annotation MQ --annotation MQRankSum --annotation ReadPosRankSum"

## Recalibrate the SNP variant quallity scores first
TYPE="SNP"


## Run the VQSR for SNPs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR VarCal for SNPs. " ' INT TERM EXIT 
${SENTIEON}/bin/sentieon driver -r ${REF} --algo VarCal -v ${SAMPLEVCF} ${RESOURCE_SNPS} ${ANNOTATE_TEXT} --var_type ${TYPE} --plot_file ${SAMPLE}.${TYPE}.plotfile --tranches_file ${SAMPLE}.${TYPE}.tranches ${SAMPLE}.${TYPE}.recal >> ${SAMPLE}.vqsr_sentieon.log 2>&1 
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi


## Apply VQSR for SNPs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR ApplyVarCal for SNPs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -r ${REF} --algo ApplyVarCal -v ${SAMPLEVCF} --var_type ${TYPE} --tranches_file ${SAMPLE}.${TYPE}.tranches --recal ${SAMPLE}.${TYPE}.recal ${SAMPLE}.${TYPE}.recaled.vcf >> ${SAMPLE}.vqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi


## Plot the report for SNP VQSR
trap 'logError " $0 stopped at line ${LINENO} Error in plot VQSR for SNPs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon plot vqsr -o ${SAMPLE}.${TYPE}.VQSR.pdf ${SAMPLE}.${TYPE}.plotfile >> ${SAMPLE}.vqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi




## Now recalibrate the INDEL variant quality scores
TYPE="INDEL"


## Run the VQSR for INDELs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR VarCal for INDELs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -r ${REF} --algo VarCal -v ${SAMPLE}.SNP.recaled.vcf ${RESOURCE_INDELS} ${ANNOTATE_TEXT} --var_type ${TYPE} --plot_file ${SAMPLE}.${TYPE}.plotfile --tranches_file ${SAMPLE}.${TYPE}.tranches ${SAMPLE}.${TYPE}.recal >> ${SAMPLE}.vqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi


## Apply VQSR for INDELs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR ApplyVarCal for INDELs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -r ${REF} --algo ApplyVarCal -v ${SAMPLE}.SNP.recaled.vcf --var_type ${TYPE} --tranches_file ${SAMPLE}.${TYPE}.tranches --recal ${SAMPLE}.${TYPE}.recal ${SAMPLE}.${TYPE}.SNP.recaled.vcf >> ${SAMPLE}.vqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi


## Plot the report for INDEL VQSR
trap 'logError " $0 stopped at line ${LINENO} Error in plot VQSR for INDELs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon plot vqsr -o ${SAMPLE}.${TYPE}.VQSR.pdf ${SAMPLE}.${TYPE}.plotfile >> ${SAMPLE}.vqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi

logInfo "[vqsr] Finished running successfully for ${SAMPLE}"
#---------------------------------------------------------------------------------------------------------------------------------------------------









#-------------------------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------------------------

## Check for the creation of final recalibrated VCF with the VQSR applied to SNPs and INDELs.
if [[ ! -s ${SAMPLE}.${TYPE}.SNP.recaled.vcf ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Output recalibrated SNP/INDEL VCF is empty."
fi

chmod g+r ${SAMPLE}.${TYPE}.SNP.recaled.vcf

#--------------------------------------------------------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------
## END
#--------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------
exit 0;

