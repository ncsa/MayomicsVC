#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## mutect.sh MANIFEST, USAGE DOCS, SET CHECKS
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

#############################################################################
#
# Script to run MuTect2 on a normal and tumor bam file
# 
#############################################################################
 USAGE:
 mutect.sh        
		   -s		<sample_name>
		   -B           <normal_bam> 
                   -T		<tumor_bam>
		   -g		<reference_genome_fasta>
                   -v           <outputVCF>
                   -I           <GATK_jar_path>
                   -J		<Java_path>
                   -M           <BCFTools_path>
		   -t		<threads>
		   -F		<shared_functions>
	           -o		<additonal options>
		   -d		Turn on debug mode
                   -h           Display this usage/help text(No arg)
                   
EXAMPLES:
mutect.sh -h
mutect.sh -s sample_name -B normal.bam -T tumor.bam -g reference_genome.fa -v output.vcf -I path/to/GATK.jar -J /path/to/java -M /path/to/BCFTools -F path/to/shared_functions -o option

NOTES: 

#############################################################################


DOCS

set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=mutect.sh
SGE_JOB_ID=TBD  # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------






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

#-------------------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------------------
## GETOPTS ARGUMENT PARSER
#------------------------------------------------------------------------------------------------------------------------------
## Check if no arguments were passed
if (($# == 0))
then
	echo -e "\nNo arguments passed.\n\n${DOCS}\n"
	exit 1
fi

## Input and Output parameters
while getopts ":hs:B:T:g:v:I:J:M:c:t:F:o:d" OPT
do
        case ${OPT} in
                h )  # Flag to dispay help message
                        echo -e "\n${DOCS}\n"
                        exit 0
			;;
	        s )  # Sample name
                        SAMPLE=${OPTARG}
                        checkArg
                        ;;
                B )  # Normal sample BAM
                        NORMAL=${OPTARG}
			checkArg
                        ;;
	        T )  # Tumor sample BAM
                        TUMOR=${OPTARG}
                        checkArg
                        ;;
                g )  # Full path to reference genome fasta file
                        REFGEN=${OPTARG}
                        checkArg
                        ;;
		v )  # Output VCF 
			OUTVCF=${OPTARG}
			checkArg
			;;
		I )  # GATK jar path
			INSTALL=${OPTARG}
			checkArg
			;;
		J )  # Java path
			JAVA=${OPTARG}
			checkArg
			;;
		M )  # BCF Tools path
			BCF=${OPTARG}
			checkArg
			;;
	        c )  # Python configuration file
                        CONFIG=${OPTARG}
                        checkArg
                        ;;
	        t )  # Number of threads
                        THR=${OPTARG}
                        checkArg
			;;
		F )  # Shared functions 
                        SHARED_FUNCTIONS=${OPTARG}
                        checkArg
                        ;;
		o )  # Extra options
                        OPTIONS=${OPTARG}
                        checkArg
                        ;;
                d )  # Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d
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

## Send Manifest to log
ERRLOG=${SAMPLE}.mutect.${SGE_JOB_ID}.log
TOOL_LOG=${SAMPLE}.mutect_tool.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 "${TOOL_LOG}"

echo "${MANIFEST}" >> "${ERRLOG}"

#SHARED_FUNCTIONS_PATH=(look at LOG_PATH in alignment)
source "${SHARED_FUNCTIONS}"

## Check if sample name is set
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Check if input files, directories, and variables are non-zero
checkVar "${NORMAL+x}" "Missing normal BAM option: -B" $LINENO
checkFile ${NORMAL} "Input normal BAM file ${NORMAL} is empty or does not exist." $LINENO

checkVar "${TUMOR+x}" "Missing tumor BAM option: -T" $LINENO
checkFile ${TUMOR} "Input tumor BAM file ${TUMOR} is empty or does not exist." $LINENO

checkVar "${REFGEN+x}" "Missing reference genome option: -g" $LINENO
checkFile ${REFGEN} "Input tumor BAM file ${REFGEN} is empty or does not exist." $LINENO

checkVar "${OUTVCF+x}" "Missing output VCF option: -v" $LINENO
#checkFile ${OUTVCF} "Output VCF file ${OUTVCF} is empty or does not exist." $LINENO

checkVar "${INSTALL+x}" "Missing GATK jar file path option: -I" $LINENO

checkVar "${JAVA+x}" "Missing Java directory option: -J" $LINENO
checkDir ${JAVA} "Reason= Java directory ${JAVA} is not a directory or does not exist." $LINENO

checkVar "${BCF+x}" "Missing BCFTools directory option: -M" $LINENO
checkDir ${BCF} "Reason= BCFTools directory ${BCF} is not a directory or does not exist." $LINENO

checkVar "${THR+x}" "Missing number of threads option: -t" $LINENO

checkVar "${SHARED_FUNCTIONS+x}" "Missing shared functions option: -F" $LINENO
checkFile ${SHARED_FUNCTIONS} "Shared functions file ${SHARED_FUNCTIONS} is empty or does not exist." $LINENO

checkVar "${OPTIONS}" "Missing additional options option: -o" $LINENO

#--------------------------------------------------------------------------------------------------------------------------------------------------

## Extra options
MUTECT_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${OPTIONS}`





#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform MuTect variant calling
#--------------------------------------------------------------------------------------------------------------------------------------------------
## Record start time
logInfo "[MuTect] START."

## first configure the MuTect run
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. MuTect2 error. Check tool log ${TOOL_LOG}. " ' INT TERM EXIT
java -jar /usr/local/apps/bioapps/gatk/gatk-3.7.0/GenomeAnalysisTK.jar \
	-T MuTect2 \
	-R ${REFGEN} \
	-I ${TUMOR} \
	-tumor tumor_${SAMPLE} \
	-I ${NORMAL} \
	-normal normal_${SAMPLE} \
	-O ${OUTVCF}
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO

checkFile ${OUTVCF} "Output somatic variant file failed to create." $LINENO

#----------------------------------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------------------------------------
## Post-Processing
#----------------------------------------------------------------------------------------------------------------------------------------------

## Clean up strelka indel output
#zcat ./strelka/results/variants/somatic.indels.vcf.gz | perl ${FIX_INDEL_GT} | gzip somatic.indels.fixed.vcf.gz
#
## Check exitCode
#

## Clean up strelka snv output
#zcat ./strelka/results/variants/somatic.snvs.vcf.gz | perl ${FIX_SNV_GT} | gzip somatic.snvs.fixed.vcf.gz
#
## Check exitCode 
#

## Combine vcfs into single output
#${BCF}/bcftools merge -m any -f PASS,. --force-sample riants/*.vcf.gz > ${SAMPLE}.vcf
#gzip ${SAMPLE}.vcf
#Not sure what this does -> | ${BCF}/bcftools plugin fill-AN-AC | ${BCF}/bcftools filter -i 'SUM(AC)>1' > ${SAMPLE}.vcf.gz

#----------------------------------------------------------------------------------------------------------------------------------------------

logInfo "Strelka workflow completed."

#----------------------------------------------------------------------------------------------------------------------------------------------
