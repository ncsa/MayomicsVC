#!/usr/bin/env bash

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
         -b 	<sorted.deduped.bam>
         -G 	<reference_genome>
         -k 	<comma,seperated,list,of,paths,to,known_sites> (omni.vcf, hapmap.vcf, indels.vcf, dbSNP.vcf)
         -I     <genomic_intervals>
         -S 	</path/to/gatk/executable> 
         -o     <extra_ApplyBQSR_options>
         -J     </path/to/java8_executable>
         -e     <java_vm_options>
         -F     </path/to/shared_functions.sh>
         -d     turn on debug mode	

 EXAMPLES:
 bqsr.sh -h
 bqsr.sh -s sample -S /path/to/gatk/executable -G reference.fa -b sorted.deduped.bam -k known1.vcf,known2.vcf,...knownN.vcf -J /path/to/java8_executable -e "'-Xms2G -Xmx8G'" -F /path/to/shared_functions.sh -o "'--createOutputBamMD5  --useOriginalQualities'" -I chr20 -d 

 NOTE: In order for getops to read in a string arguments for -o (extra_ApplyBQSR_options) or -e (java_vm_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.

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


while getopts ":hs:S:G:b:k:J:e:F:o:I:d" OPT
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
		S ) # Full path to gatk executable 
			GATKEXE=${OPTARG}
			checkArg
			;;
		G ) # Full path to reference fasta
			REF=${OPTARG}
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
        J )  # Path to JAVA8 exectable. The variable needs to be small letters so as not to explicitly change the user's $PATH variable 
             java=${OPTARG}
             checkArg
             ;;
        e )  # JAVA options string to pass into the gatk command 
             JAVA_OPTS_STRING=${OPTARG}
             checkArg
             ;;
		F )  # Path to shared_functions.sh
		     SHARED_FUNCTIONS=${OPTARG}
			 checkArg
			;;
        o ) # Extra options and arguments to ApplyBQSR, input as a long string, can be empty if desired
              APPLYBQSR_OPTIONS=${OPTARG}
              checkArg
              ;;
        I )  # Genomic intervals overwhich to operate
              INTERVALS=${OPTARG}
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
truncate -s 0 ${SAMPLE}.bqsr_gatk.log


## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check java path and options 
checkVar "${java+x}" "Missing JAVA path option: -J" $LINENO
checkFileExe ${java} "REASON=JAVA file ${java} is not executable or does not exist." $LINENO
checkVar "${JAVA_OPTS_STRING+x}" "Missing specification of JAVA memory options: -e" $LINENO

## Check if the GATK executable option was passed in.
checkVar "${GATKEXE+x}" "Missing GATK path option: -S" $LINENO

## Check if the GATK executable is present.
checkFileExe ${GATKEXE} "REASON=GATK file ${GATKEXE} is not executable or does not exist." $LINENO 

## Check if the reference option was passed in
checkVar "${REF+x}" "Missing reference genome option: -G" $LINENO

## Check if the reference fasta file is present.
checkFile ${REF} "Reference genome file ${REF} is empty or does not exist." $LINENO

## Check if the BAM input file option was passed in
checkVar "${INPUTBAM+x}" "REASON=Missing input BAM option: -b" $LINENO

## Check if the BAM input file is present.
checkFile ${INPUTBAM} "Input BAM ${INPUTBAM} is empty or does not exist." $LINENO
checkFile ${INPUTBAM}.bai "Input BAM index ${INPUTBAM}.bai is empty or does not exist." $LINENO

## Check if the known sites file option is present.
checkVar "${KNOWN+x}" "Missing known sites option: -k" $LINENO

checkVar "${APPLYBQSR_OPTIONS+x}" "Missing extra ApplyBQSR options option: -o" $LINENO

checkVar "${INTERVALS+x}" "Missing Intervals option: -I" $LINENO

#---------------------------------------------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------------------------------------------
## FILENAME AND OPTION PARSING
#---------------------------------------------------------------------------------------------------------------------------

## Parse known sites list of multiple files. Create multiple -k flags for gatk
KNOWNSITES=$( echo --known-sites ${KNOWN} | sed "s/,/ --known-sites /g" | tr "\n" " " )

## Extra options
APPLYBQSR_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${APPLYBQSR_OPTIONS}`
JAVA_OPTS_PARSED=`sed -e "s/'//g" <<< ${JAVA_OPTS_STRING}`

## Defining file names
TOOL_LOG=${SAMPLE}.bqsr_gatk.log


#---------------------------------------------------------------------------------------------------------------------------
## Perform bqsr with GATK 
#---------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[bqsr] START. Generating the bqsr model"


#Calculate required modification of the quality scores in the BAM
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in bqsr Step1: Calculate required modification of the quality scores in the BAM. " ' INT TERM EXIT
${GATKEXE} --java-options "${JAVA_OPTS_PARSED}" BaseRecalibrator --reference ${REF} --input ${INPUTBAM} ${KNOWNSITES} --output ${SAMPLE}.${INTERVALS}.recal_data.table --intervals ${INTERVALS} >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO
logInfo "[bqsr] Finished generating the bqsr table for ${SAMPLE}.${INTERVALS}" 

## Record start time
logInfo "[bqsr] START. Generate the bqsr'd bam file"


#Calculate required modification of the quality scores in the BAM
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in bqsr Step2: Generate a BAM with modifications of the quality scores. " ' INT TERM EXIT
${GATKEXE} --java-options "${JAVA_OPTS_PARSED}" ApplyBQSR --reference ${REF} --input ${INPUTBAM} --output ${SAMPLE}.${INTERVALS}.bam -bqsr ${SAMPLE}.${INTERVALS}.recal_data.table ${APPLYBQSR_OPTIONS_PARSED} --intervals ${INTERVALS} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30  >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO
logInfo "[bqsr] Finished running successfully and generated the bam file ${SAMPLE}.${INTERVALS}.bam" 





#---------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#---------------------------------------------------------------------------------------------------------------------------

# Check for the creation of the bam file for input to Haplotyper. Open read permissions for the group.

checkFile ${SAMPLE}.${INTERVALS}.bam "Recalibrated file ${SAMPLE}.${INTERVALS}.bam is empty." $LINENO
checkFile ${SAMPLE}.${INTERVALS}.bai "Output recalibrated BAM index ${SAMPLE}.${INTERVALS}.bai is empty." ${LINENO}
chmod g+r ${SAMPLE}.${INTERVALS}.bam
chmod g+r ${SAMPLE}.${INTERVALS}.bai

#---------------------------------------------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
## END
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
exit 0;
