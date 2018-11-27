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
	 -G 	<reference_genome>
	 -t 	<threads>
	 -b 	<sorted.deduped.realigned.bam>
	 -k 	<known_sites> (omni.vcf, hapmap.vcf, indels.vcf, dbSNP.vcf)
         -e     </path/to/env_profile_file>
	 -F     </path/to/shared_functions.sh>
	 -d	turn on debug mode	

 EXAMPLES:
 bqsr.sh -h
 bqsr.sh -s sample -S /path/to/sentieon_directory -G reference.fa -t 12 -b sorted.deduped.realigned.bam -k known1.vcf,known2.vcf,...knownN.vcf -e /path/to/env_profile_file -F /path/to/shared_functions.sh -d 

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


while getopts ":hs:S:G:t:b:k:e:F:d" OPT
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
		S ) # Full path to Sentieon
			SENTIEON=${OPTARG}
			checkArg
			;;
		G ) # Full path to reference fasta
			REF=${OPTARG}
			checkArg
			;;
		t ) # Number of threads available
			NTHREADS=${OPTARG}
			checkArg
			;;
		b ) # Full path to DeDuped BAM used as input
			INPUTBAM=${OPTARG}
			checkArg
			;;
		k ) # Full path to known sites files
			KNOWN=${OPTARG}
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
		d ) # Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d.
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

source ${SHARED_FUNCTIONS}

## Check if sample name is present.
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.bqsr.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.bqsr_sentieon.log


## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with environmental profile variables
checkVar "${ENV_PROFILE+x}" "Missing environmental profile option: -e" $LINENO
source ${ENV_PROFILE}

## Check if the Sentieon executable option was passed in.
checkVar "${SENTIEON+x}" "Missing Sentieon path option: -S" $LINENO

## Check if the Sentieon executable is present.
checkDir ${SENTIEON} "Sentieon directory ${SENTIEON} is not a directory or does not exist." $LINENO 

#Check if the number of threads is present.
checkVar "${NTHREADS+x}" "Missing threads option: -t" $LINENO

## Check if the reference option was passed in
checkVar "${REF+x}" "Missing reference genome option: -G" $LINENO

## Check if the reference fasta file is present.
checkFile ${REF} "Reference genome file ${REF} is empty or does not exist." $LINENO

## Check if the BAM input file option was passed in
checkVar "${INPUTBAM+x}" "REASON=Missing input BAM option: -b" $LINENO

## Check if the BAM input file is present.
checkFile ${INPUTBAM} "Input BAM ${INPUTBAM} is empty or does not exist." $LINENO
checkFile ${INPUTBAM}.bai "Input BAM index ${INPUTBAM} is empty or does not exist." $LINENO

## Check if the known sites file option is present.
checkVar "${KNOWN+x}" "Missing known sites option ${KNOWN}: -k" $LINENO

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


#Calculate required modification of the quality scores in the BAM
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in bqsr Step1: Calculate required modification of the quality scores in the BAM. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} --algo QualCal -k ${SPLITKNOWN} ${SAMPLE}.recal_data.table >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

#Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in bqsr Step2: Apply the recalibration to calculate the post calibration data table and additionally apply the recalibration on the BAM file. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} -q ${SAMPLE}.recal_data.table --algo QualCal -k ${SPLITKNOWN} ${SAMPLE}.recal_data.table.post >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

#Create data for plotting
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in bqsr Step3: Create data for plotting. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} --algo QualCal --plot --before ${SAMPLE}.recal_data.table --after ${SAMPLE}.recal_data.table.post ${SAMPLE}.recal.csv >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

#Plot the calibration data tables, both pre and post, into graphs in a pdf
TRAP_LINE=$(($LINENO + 1))
trap 'logError "$0 stopped at line ${TRAP_LINE}. Error in bqsr Step4: Plot the calibration data tables, both pre and post, into graphs in a pdf. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon plot bqsr -o ${SAMPLE}.recal_plots.pdf ${SAMPLE}.recal.csv >> ${SAMPLE}.bqsr_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

logInfo "[bqsr] Finished running successfully for ${SAMPLE}" 
#------------------------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#------------------------------------------------------------------------------------------------------------------------------------

# Check for the creation of the recal_data.table necessary for input to Haplotyper. Open read permissions for the group.
# The other files created in BQSR are not necessary for the workflow to run, so I am not performing checks on them.

checkFile ${SAMPLE}.recal_data.table "Recal data table ${SAMPLE}.recal_data.table is empty." $LINENO

chmod g+r ${SAMPLE}.recal_data.table

#-----------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;

