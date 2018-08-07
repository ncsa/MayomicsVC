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
	 -G 	<reference_genome>
	 -t 	<threads>
	 -b 	<sorted.deduped.realigned.bam>
	 -k 	<known_sites> (omni.vcf, hapmap.vcf, indels.vcf, dbSNP.vcf)
	 -d	turn on debug mode	

 EXAMPLES:
 bqsr.sh -h
 bqsr.sh -s sample -S /path/to/sentieon_directory -L sentieon_license_number -G reference.fa -t 12 -b sorted.deduped.realigned.bam -k known1.vcf,known2.vcf,...knownN.vcf -d 

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

function checkArg()
{
    if [[ "${OPTARG}" == -* ]]; then
        echo -e "\nError with option -${OPT} in command. Option passed incorrectly or without argument.\n"
        echo -e "\n${DOCS}\n"
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


while getopts ":hs:S:L:G:t:b:k:d" OPT
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
		G ) # Full path to reference fasta. String variable invoked with -G
			REF=${OPTARG}
			checkArg
			;;
		t ) # Number of threads available. Integer invoked with -t
			NTHREADS=${OPTARG}
			checkArg
			;;
		b ) # Full path to DeDuped BAM used as input. String variable invoked with -b
			INPUTBAM=${OPTARG}
			checkArg
			;;
		k ) # Full path to known sites files. String variable invoked with -k
			KNOWN=${OPTARG}
			checkArg
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
if [[ -z ${SAMPLE+x} ]] ## NOTE: ${VAR+x} is used for variable expansions, preventing unset variable error from set -o nounset. When $VAR is not set, we set it to "x" and throw the error.
then
	echo -e "$0 stopped at line ${LINENO}. \nREASON=Missing sample name option: -s"
	exit 1
fi

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.bqsr.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.bqsr_sentieon.log


## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if the Sentieon executable option was passed in.
if [[ -z ${SENTIEON+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon path option: -S"
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

## Check if the reference option was passed in
if [[ -z ${REF+x} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Missing reference genome option: -G"
fi

## Check if the reference fasta file is present.
if [[ ! -s ${REF} ]]
then
	EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REF} is empty or does not exist."
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
	logError "$0 stopped at line $LINENO. \nREASON=Input BAM ${INPUTBAM} is empty or does not exist."
fi

if [[ ! -s ${INPUTBAM}.bai ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Input BAM index ${INPUTBAM} is empty or does not exist."
fi


## Check if the known sites file option is present.
if [[ -z ${KNOWN+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing known sites option ${KNOWN}: -k"
fi

## Check if Sentieon license string is present.
if [[ -z ${LICENSE+x} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Missing Sentieon liscence required option: -L"
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------------------------------------------------------------------
## FILENAME AND OPTION PARSING
#---------------------------------------------------------------------------------------------------------------------------------------------------

## Parse known sites list of multiple files. Create multiple -k flags for sentieon
SPLITKNOWN=`sed -e 's/,/ -k /g' <<< ${KNOWN}`

#---------------------------------------------------------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform bqsr with Sention
#--------------------------------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[bqsr] START. Performing bqsr on the input BAM to produce bqsr table."

export SENTIEON_LICENSE=${LICENSE}

#Calculate required modification of the quality scores in the BAM

trap 'logError " $0 stopped at line ${LINENO}. Error in bqsr Step1: Calculate required modification of the quality scores in the BAM. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} --algo QualCal -k ${SPLITKNOWN} ${SAMPLE}.recal_data.table >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi


#Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file
trap 'logError " $0 stopped at line ${LINENO}. Error in bqsr Step2: Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${SAMPLE}.recal_data.table --algo QualCal -k ${SPLITKNOWN} ${SAMPLE}.recal_data.table.post >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}."
fi	


#Create data for plotting
trap 'logError " $0 stopped at line ${LINENO}. Error in bqsr Step3: Create data for plotting. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} --algo QualCal --plot --before ${SAMPLE}.recal_data.table --after ${SAMPLE}.recal_data.table.post ${SAMPLE}.recal.csv >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in bqsr Step3: Create data for plotting"
fi	


#Plot the calibration data tables, both pre and post, into graphs in a pdf
trap 'logError "$0 stopped at line ${LINENO}. Error in bqsr Step4: Plot the calibration data tables, both pre and post, into graphs in a pdf. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon plot bqsr -o ${SAMPLE}.recal_plots.pdf ${SAMPLE}.recal.csv >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in bqsr Step4: Plot the calibration data tables, both pre and post, into graphs in a pdf"
fi	

logInfo "[bqsr] Finished running successfully for ${SAMPLE}" 
#------------------------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#------------------------------------------------------------------------------------------------------------------------------------

# Check for the creation of the recal_data.table necessary for input to Haplotyper. Open read permissions for the group.
# The other files created in BQSR are not necessary for the workflow to run, so I am not performing checks on them.

if [[ ! -s ${SAMPLE}.recal_data.table ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Recal data table ${SAMPLE}.recal_data.table is empty."
fi

chmod g+r ${SAMPLE}.recal_data.table

#-----------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;

