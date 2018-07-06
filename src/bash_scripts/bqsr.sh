#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------------
## bqsr.sh MANIFEST, USAGE DOCS, SET CHECKS
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

############################################################################################################################
#
# Perform Base Quality Score Recalibration (bqsr) on the bam produced in the Deduplication stage of the Mayomics workflow.
# Step 1/3 in Single Sample Variant Calling.
#
############################################################################################################################

 USAGE:
 bqsr.sh -s 	<sample_name>
	 -S 	</path/to/sentieon> 
	 -L	<sentieon_license>
	 -G 	</path/to/ref.fa>
	 -t 	<threads>
	 -b 	</path/to/Sorted_deDuped.bam>
	 -D 	</path/to/dbsnp.vcf>
	 -k 	</path/to/known_indels.vcf>
	 -d	turn on debug mode	

 EXAMPLES:
 bqsr.sh -h
 bqsr.sh -s sample -S sentieon -L sentieon_License -G ref.fa -t 12 -b sample.bam -D dbsnp.vcf -k known_indels.vcf -d 

############################################################################################################################

DOCS

set -o errexit
set -o pipefail
set -o nounset


SCRIPT_NAME=bqsr.sh
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
  
    if [[ -n ${ERRLOG+x} ]]; then
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

    if [[ -z ${EXITCODE+x} ]]
    then
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


while getopts ":hs:S:L:G:t:b:D:k:d" OPT
do
	case ${OPT} in
		h ) # flag to display help message
			echo -e "\n${DOCS}\n "
			exit 0; 
			;;
		s ) # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			#echo ${SAMPLE}
			;;
		S ) # Full path to Sentieon executable. String variable invoked with -S
			SENTIEON=${OPTARG}
			#echo ${SENTIEON}
			;;
		L ) # Sentieon license number. Invoked with -L 
			LICENSE=${OPTARG}
			#echo ${LICENSE}
			;;
		G ) # Full path to reference fasta. String variable invoked with -r
			REF=${OPTARG}
			#echo ${REF}
			;;
		t ) # Number of threads available. Integer invoked with -t
			NTHREADS=${OPTARG}
			#echo ${NTHREADS}
			;;
		b ) # Full path to DeDuped BAM used as input. String variable invoked with -i
			INPUTBAM=${OPTARG}
			#echo ${INPUTBAM}
			;;
		D ) # Full path to DBSNP file. String variable invoked with -D.
			DBSNP=${OPTARG}
			#echo ${DBSNP}
			;;
		k ) # Full path to known site indel file (dbSNP and known indels VCF), separated by a comma no space. String variable invoked with -k #MIGHT TAKE IN MULTILPE FILES
			KNOWN=${OPTARG}
			#echo ${KNOWN}
			;;
		d ) # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text.
			echo -e "\nDebug mode is ON.\n"
			set -x
			;;
		\? ) # Check for unsupported flag, print usage and exit
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
## Check if sample name is present.
if [[ -z ${SAMPLE+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=String for sample name is not present."
fi

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.bqsr.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"

## Create log for the specific tool in this script, bqsr. 
TOOL_LOG=${SAMPLE}.bqsr_sentieon.log
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
if [[ ! -s ${SENTIEON} ]]
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

## Check if the BAM input file option was passed in
if [[ -z ${INPUTBAM+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing input BAM required option: -b"
fi

## Check if the BAM input file is present.
if [[ ! -s ${INPUTBAM} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Input BAM ${INPUTBAM} is not present or does not exist."
fi

## Check if dbSNP input file option was passed in
if [[ -z ${DBSNP+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing dbSNP required option: -D"
fi

## Check if dbSNP file is present.
if [[ ! -f ${DBSNP} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=DBSNP ${DBSNP} is not present or does not exist."
fi

## Check if the known indels file is present.
if [[ -z ${KNOWN+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing known indels file required option: -k"
fi

## Check if the known indels file is present.
if [[ ! -f ${KNOWN} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Known indels file ${KNOWN} is not present or does not exist."
fi

## Check if Sentieon license string is present.
if [[ -z ${LICENSE+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon liscence required option: -L"
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform bqsr with Sention
#--------------------------------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[bqsr] START. Performing bqsr on the input BAM to produce bqsr table."

export SENTIEON_LICENSE=${LICENSE}

#Calculate required modification of the quality scores in the BAM

trap 'logError " $0 stopped at line ${LINENO}. Error in bqsr Step1: Calculate required modification of the quality scores in the BAM. " ' INT TERM EXIT
${SENTIEON} driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} --algo QualCal -k ${DBSNP} -k ${KNOWN} ${SAMPLE}.recal_data.table >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
	exit ${EXITCODE};
fi


#Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file
trap 'logError " $0 stopped at line ${LINENO}. Error in bqsr Step2: Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file. " ' INT TERM EXIT
${SENTIEON} driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${SAMPLE}.recal_data.table --algo QualCal -k ${DBSNP} -k ${KNOWN} ${SAMPLE}.recal_data.table.post >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
	exit ${EXITCODE};
fi	


#Create data for plotting
trap 'logError " $0 stopped at line ${LINENO}. Error in bqsr Step3: Create data for plotting. " ' INT TERM EXIT
${SENTIEON} driver -t ${NTHREADS} --algo QualCal --plot --before ${SAMPLE}.recal_data.table --after ${SAMPLE}.recal_data.table.post ${SAMPLE}.recal.csv >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in bqsr Step3: Create data for plotting"
	exit ${EXITCODE};
fi	


#Plot the calibration data tables, both pre and post, into graphs in a pdf
trap 'logError "$0 stopped at line ${LINENO}. Error in bqsr Step4: Plot the calibration data tables, both pre and post, into graphs in a pdf. " ' INT TERM EXIT
${SENTIEON} plot bqsr -o ${SAMPLE}.recal_plots.pdf ${SAMPLE}.recal.csv >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in bqsr Step4: Plot the calibration data tables, both pre and post, into graphs in a pdf"
	exit ${EXITCODE};
fi	
#------------------------------------------------------------------------------------------------------------------------------------



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

