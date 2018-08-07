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
	 -H <hapmap.vcf>
	 -O <omni.vcf>
	 -T <1000genomes.vcf>
	 -D <dbsnp.vcf>
	 -m <mills.vcf>
	 -d turn on debug mode

 EXAMPLES:
 vqsr.sh -h
 vqsr.sh -s sample -S sentieon -L sentieon_license -G ref.fa -V sample.vcf -H hapmap.vcf -O omni.vcf -T thousandg.vcf -D dbsnp.vcf -m mills.vcf -d

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
while getopts ":hs:S:L:G:V:H:O:T:D:m:d" OPT
do
	case ${OPT} in
		h ) # Flag to display usage
			echo -e "\n${DOCS}\n"
			exit 0
			;;
		s ) # Sample name. Invoked with -s.
			SAMPLE=${OPTARG}
			checkArg
			;;
		S ) #Full path to Sentieon. Invoked with -S.
			SENTIEON=${OPTARG}
			checkArg
			;;
		L ) #Sentieon license. Invoked with -L.
			LICENSE=${OPTARG}
			checkArg
			;;
		G ) #Reference genome. Invoked with -G.
			REF=${OPTARG}
			checkArg
			;;
		V ) #Sample VCF file output from Haplotyper. Invoked with -V.
			SAMPLEVCF=${OPTARG}
			checkArg
			;;
		H ) #Hapmap VCF file as a known site. Invoked with -H
			HAPMAP=${OPTARG}
			checkArg
			;;
		O ) #Omni VCF file as a known site. Invoked with -O.
			OMNI=${OPTARG}
			checkArg
			;;
		T ) #1000 genomes VCF file as a known site. Invoked with -T.
			THOUSANDG=${OPTARG}
			checkArg
			;;
		D ) #dbSNP VCF file as a known site. Invoked with -D.
			DBSNP=${OPTARG}
			checkArg
			;;
		m ) #Mills VCF file as a known site. Invoked with -m.
			MILLS=${OPTARG}
			checkArg
			;;
		d ) # Turn on debug mode. Initiates 'set -x' to print all text.
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

## Create log for the specific tool in this script, bqsr. 
TOOL_LOG=${SAMPLE}.vqsr_sentieon.log
truncate -s 0 "${TOOL_LOG}"

## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if the Sentieon executable option was passed in.
if [[ -z ${SENTIEON+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon executable required option: -S"
fi

## Check if the Sentieon executable is present.
if [[ ! -d ${SENTIEON} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Sentieon executable ${SENTIEON} is not present or does not exist."
fi

## Check if the Sentieon license option was passed in.
if [[ -z ${LICENSE+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon license required option: -S"
fi

## Check if the reference option was passed in
if [[ -z ${REF+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing reference fasta file required option: -G"
fi

## Check if the reference fasta file is present.
if [[ ! -s ${REF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome fasta file ${REF} is not present or does not exist."
fi

## Check if the sample VCF input file option was passed in.
if [[ -z ${SAMPLEVCF+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing sample input VCF required option: -V"
fi

## Check if the sample VCF input file is present.
if [[ ! -s ${SAMPLEVCF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input sample vcf ${SAMPLEVCF} is not present or does not exist."
fi

## Check if the Hapmap VCF input file option was passed in.
if [[ -z ${HAPMAP+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Hapmap input VCF required option: -H"
fi

## Check if the Hapmap VCF input file is present.
if [[ ! -s ${HAPMAP} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Hapmap VCF ${HAPMAP} is not present or does not exist."
fi

## Check if the Omni VCF input file option was passed in.
if [[ -z ${OMNI+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Omni input VCF required option: -O"
fi

## Check if the Omni VCF input file is present.
if [[ ! -s ${OMNI} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Omni VCF ${OMNI} is not present or does not exist."
fi

## Check if the 1000 genomes VCF input file option was passed in.
if [[ -z ${THOUSANDG+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing 1000 genomes input VCF required option: -T"
fi

## Check if the 1000 genomes VCF input file is present.
if [[ ! -s ${THOUSANDG} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=1000 genomes VCF ${THOUSANDG} is not present or does not exist."
fi

## Check if the dbSNP VCF input file option was passed in.
if [[ -z ${DBSNP+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing DNPSNP input VCF required option: -D"
fi

## Check if the dbSNP VCF input file is present.
if [[ ! -s ${DBSNP} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=dbSNP VCF ${DBSNP} is not present or does not exist."
fi


## Check if the Mills VCF input file option was passed in.
if [[ -z ${MILLS+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing Mills input VCF required option: -m"
fi

## Check if the Mills VCF input file is present.
if [[ ! -s ${MILLS} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Mills VCF ${MILLS} is not present or does not exist."
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

## Create the resource argument
RESOURCE_TEXT="--resource ${THOUSANDG} --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
RESOURCE_TEXT="${RESOURCE_TEXT} --resource ${OMNI} --resource_param omni,known=false,training=true,truth=false,prior=12.0 "
RESOURCE_TEXT="${RESOURCE_TEXT} --resource ${DBSNP} --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
RESOURCE_TEXT="${RESOURCE_TEXT} --resource ${HAPMAP} --resource_param hapmap,known=false,training=true,truth=true,prior=15.0"


## Run the VQSR for SNPs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR VarCal for SNPs. " ' INT TERM EXIT 
${SENTIEON}/bin/sentieon driver -r ${REF} --algo VarCal -v ${SAMPLEVCF} ${RESOURCE_TEXT} ${ANNOTATE_TEXT} --var_type ${TYPE} --plot_file ${SAMPLE}.${TYPE}.plotfile --tranches_file ${SAMPLE}.${TYPE}.tranches ${SAMPLE}.${TYPE}.recal
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi


## Apply VQSR for SNPs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR ApplyVarCal for SNPs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -r ${REF} --algo ApplyVarCal -v ${SAMPLEVCF} --var_type ${TYPE} --tranches_file ${SAMPLE}.${TYPE}.tranches --recal ${SAMPLE}.${TYPE}.recal ${SAMPLE}.${TYPE}.recaled.vcf
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi


## Plot the report for SNP VQSR
trap 'logError " $0 stopped at line ${LINENO} Error in plot VQSR for SNPs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon plot vqsr -o ${SAMPLE}.${TYPE}.VQSR.pdf ${SAMPLE}.${TYPE}.plotfile
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi




## Now recalibrate the INDEL variant quality scores
TYPE="INDEL"

## Create the resource argument
RESOURCE_TEXT=""
RESOURCE_TEXT="--resource ${MILLS} --resource_param Mills,known=false,training=true,truth=true,prior=12.0 "
RESOURCE_TEXT="${RESOURCE_TEXT} --resource ${DBSNP} --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0"

## Run the VQSR for INDELs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR VarCal for INDELs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -r ${REF} --algo VarCal -v ${SAMPLE}.SNP.recaled.vcf ${RESOURCE_TEXT} ${ANNOTATE_TEXT} --var_type ${TYPE} --plot_file ${SAMPLE}.${TYPE}.plotfile --tranches_file ${SAMPLE}.${TYPE}.tranches ${SAMPLE}.${TYPE}.recal
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi


## Apply VQSR for INDELs
trap 'logError " $0 stopped at line ${LINENO} Error in VQSR ApplyVarCal for INDELs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -r ${REF} --algo ApplyVarCal -v ${SAMPLE}.SNP.recaled.vcf --var_type ${TYPE} --tranches_file ${SAMPLE}.${TYPE}.tranches --recal ${SAMPLE}.${TYPE}.recal ${SAMPLE}.${TYPE}.SNP.recaled.vcf
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi


## Plot the report for INDEL VQSR
trap 'logError " $0 stopped at line ${LINENO} Error in plot VQSR for INDELs. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon plot vqsr -o ${SAMPLE}.${TYPE}.VQSR.pdf ${SAMPLE}.${TYPE}.plotfile
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
        logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
        exit ${EXITCODE};
fi
#---------------------------------------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#------------------------------------------------------------------------------------------------------------------------------------


logInfo "[bqsr] Finished running successfully for ${SAMPLE}"
#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;

