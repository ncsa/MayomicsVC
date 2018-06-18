#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------------
## BQSR.sh MANIFEST, USAGE DOCS, SET CHECKS
#------------------------------------------------------------------------------------------------------------------------------


read -r -d '' MANIFEST << MANIFEST
*******************************************
`readlink -m $0` was called by: `whoami` on `date`
command line input: ${@}
*******************************************
MANIFEST
echo "${MANIFEST}"

read -r -d '' DOCS << DOCS

############################################################################################################################
#
# Perform Base Quality Score Recalibration (BQSR) on the bam produiced in the Deduplication stage of the Mayomics workflow.
# Step 1/3 in Single Sample Variant Calling.
#
############################################################################################################################

 USAGE:
 BQSR.sh -s 	<sample_name>
	 -O 	</path/to/output_dir>
	 -S 	</path/to/sentieon> 
	 -L	<sentieon_license>
	 -G 	</path/to/ref.fa>
	 -t 	<threads>
	 -b 	</path/to/deDuped.bam>
	 -D 	</path/to/dbsnp.vcf>
	 -k 	</path/to/known_indels.vcf>
	 -e 	</path/to/error_log>
	 -d 	debug_mode [false]

 EXAMPLES:
 BQSR.sh -h
 BQSR -s sample -O output_dir -S sentieon -r ref.fa -t 12 -i sample.bam -D dbsnp.vcf -k known_indels.vcf -e error.log 

############################################################################################################################

DOCS

set -o errexit
set -o pipefail
#set -o nounset


SCRIPT_NAME=BQSR.sh
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
while getopts ":hs:O:S:L:G:t:b:D:k:e:d" OPT
do
	case ${OPT} in
		h ) # flag to display help message
			echo " "
			echo "Usage:"
			echo " "
			echo "  bash BQSR.sh -h       Display this help message." 
			echo "  bash BQSR.sh [-s <sample_name>] [-O </path/to/output_dir>] [-S </path/to/sentieon>] [-G </path/to/ref.fa>] [-t <threads>] [-b </path/to/deDuped.bam>] [-D </path/to/dbsnp.vcf>] [-k </path/to/known_indels.vcf>] [-e </path/to/error_log>] [-d debug_mode [false]]"
			echo " "
			exit 0; 
			;;
		s ) # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			#echo ${SAMPLE}
			;;
		O ) # Output directory. String variable invoked with -O
			OUTDIR=${OPTARG}
			#echo ${OUTDIR}
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
		e ) # Full path to error log file. String variable invoked with -e
			ERRLOG=${OPTARG}
			#echo ${ERRLOG}
			;;
		d ) # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text.
			DEBUG=${OPTARG}
			#echo ${DEBUG}
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

## Check if the known indels file is present.
if [[ ! -f ${KNOWN} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Known indels file ${KNOWN} is not present or does not exist."
fi

## Check if Sentieon license string is present.
if [[ -z ${LICENSE} ]]
then
	EXITCODE=1
	logError "$0 stopped at line $LINENO. \nREASON=Sentieon license ${LICENSE} is not present or does not exist."
fi
#--------------------------------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform BQSR with Sention
#--------------------------------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[BQSR] START."

export SENTIEON_LICENSE=${LICENSE}

#Calculate required modification of the quality scores in the BAM
${SENTIEON} driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} --algo QualCal -k ${DBSNP} -k ${KNOWN} ${OUTDIR}/${SAMPLE}.recal_data.table
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in BQSR Step1: Calculate required modification of the quality scores in the BAM."
	exit ${EXITCODE};
fi


#Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file
${SENTIEON} driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${OUTDIR}/${SAMPLE}.recal_data.table --algo QualCal -k ${DBSNP} -k ${KNOWN} ${OUTDIR}/${SAMPLE}.recal_data.table.post	
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in BQSR Step2: Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file."
	exit ${EXITCODE};
fi	


#Create data for plotting
${SENTIEON} driver -t ${NTHREADS} --algo QualCal --plot --before ${OUTDIR}/${SAMPLE}.recal_data.table --after ${OUTDIR}/${SAMPLE}.recal_data.table.post ${OUTDIR}/${SAMPLE}.recal.csv
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in BQSR Step3: Create data for plotting"
	exit ${EXITCODE};
fi	


#Plot the calibration data tables, both pre and post, into graphs in a pdf
${SENTIEON} plot bqsr -t ${NTHREADS} -o ${OUTDIR}/${SAMPLE}.recal_plots.pdf ${OUTDIR}/${SAMPLE}.recal.csv
EXITCODE=$?
if [[ ${EXITCODE} -ne 0 ]]
then
	logError "$0 stopped at line ${LINENO} with exit code ${EXITCODE}. Error in BQSR Step4: Plot the calibration data tables, both pre and post, into graphs in a pdf"
	exit ${EXITCODE};
fi	
#------------------------------------------------------------------------------------------------------------------------------------

## Open read permissions to the user group
chmod g+r -R ${OUTDIR}


logInfo "[BQSR] Finished running successfully. Output for ${SAMPLE} can be found in ${OUTDIR}/" 
#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;

