#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## combine_variants.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Script to merge the output VCFs from Strelka and MuTect.
# Part of the MayomicsVC Workflow.
# 
#############################################################################
 USAGE:
combine_variants.sh        
                       -s           <sample_name>
                       -S           <Strelka_vcf>
                       -M           <MuTect_vcf>           
                       -g           <reference_genome_fasta>
                       -G           <GATK_jar_path>
                       -J           <Java_path>
                       -Z           <bgzip_path>
                       -t           <threads>
                       -F           <shared_functions>
                       -e           <environmental_profile>
                       -p           <prioritization option string>   
                       -o           <additonal options>
                       -d           Turn on debug mode
                       -h           Display this usage/help text(No arg)
                   
EXAMPLES:
combine_variants.sh -h
combine_variants.sh -s sample_name -S strelka.vcf -M mutect.vcf -g reference_genome.fa -G /path/to/GATK -J /path/to/java -Z /path/to/bgzip -F path/to/shared_functions -e /path.to/environmentprofile.file -o "'--genotypeMergeOptions PRIORITIZE -priority mutect.vcf,strelka.vcf"

NOTES: 
Common PrioritizationOptionString would be "strelka,mutect" or "mutect,strelka". 
Is a comma-separated list of those two tools, no CAPS. 

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
while getopts ":hs:S:M:g:G:J:Z:t:F:e:p:o:d" OPT
do
	case ${OPT} in
		h )  # Flag to display help message
			echo -e "\n${DOCS}\n"
			exit 0
			;;
		s )  # Sample name
			SAMPLE=${OPTARG}
			checkArg
			;;
		S )  # Strelka output vcf
			STRELKA_VCF=${OPTARG}
			checkArg
			;;
		M )  # MuTect ouput vcf
			MUTECT_VCF=${OPTARG}
			checkArg
			;;
		g )  # Full path to reference genome fasta file
			REFGEN=${OPTARG}
			checkArg
			;;
		G )  # GATK path
			GATK=${OPTARG}
			checkArg
			;;
		J )  # Java path
			JAVA=${OPTARG}
			checkArg
			;;
		Z )  # Bgzip path
			BGZIP=${OPTARG}
			checkArg
			;;
		t )  # Number of threads
			THR=${OPTARG}
			checkArg
			;;
		F )  # Path to shared_functions.sh
			SHARED_FUNCTIONS=${OPTARG}
			checkArg
			;;
                e )  # Path to file with environmental profile variables
                        ENV_PROFILE=${OPTARG}
                        checkArg
                        ;;
                p )  # Prioritization string
                        PRIORITIZATION_STRING=${OPTARG}
                        checkArg
                        ;;
		o )  # Extra options string
			OPTIONS=${OPTARG}
			checkArg
			;;
		d )  # Turn on debug mode. Initiates 'set -x' to print all code
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

source "${SHARED_FUNCTIONS}"

## Create log for JOB_ID/script and tool
ERRLOG=${SAMPLE}.combine_variants.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=${SAMPLE}.CombineVariants.log
truncate -s 0 ${TOOL_LOG}

## Send Manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with environmental profile variables
checkVar "${ENV_PROFILE+x}" "Missing environmental profile option: -e" $LINENO                                       
source ${ENV_PROFILE}



## Check if sample name is set
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Check if input files, directories, and variables are non-zero
checkVar "${STRELKA_VCF+x}" "Missing Strelka input vcf option: -S" $LINENO
checkFile ${STRELKA_VCF} "Input Strelka vcf file ${STRELKA_VCF} is empty or does not exist." $LINENO

checkVar "${MUTECT_VCF+x}" "Missing MuTect input vcf option: -M" $LINENO
checkFile ${MUTECT_VCF} "Input MuTect vcf file ${MUTECT_VCF} is empty or does not exist." $LINENO

checkVar "${REFGEN+x}" "Missing reference genome option: -g" $LINENO
checkFile ${REFGEN} "Input tumor BAM file ${REFGEN} is empty or does not exist." $LINENO

checkVar "${GATK+x}" "Missing GATK directory path option: -G" $LINENO
checkDir ${GATK} "Reason= GATK directory ${GATK} is not a directory or does not exist." $LINENO

checkVar "${JAVA+x}" "Missing Java directory option: -J" $LINENO
checkDir ${JAVA} "Reason= Java directory ${JAVA} is not a directory or does not exist." $LINENO

checkVar "${BGZIP+x}" "Missing bgzip directory option: -Z" $LINENO
checkDir ${BGZIP} "Reason= bgzip directory ${BGZIP} is not a directory or does not exist." $LINENO

checkVar "${THR+x}" "Missing number of threads option: -t" $LINENO

checkVar "${SHARED_FUNCTIONS+x}" "Missing shared functions option: -F" $LINENO
checkFile ${SHARED_FUNCTIONS} "Shared functions file ${SHARED_FUNCTIONS} is empty or does not exist." $LINENO

checkVar "$PRIORITIZATION_STRING{}" "Missing prioritization string option: -p" $LINENO
checkVar "${OPTIONS}" "Missing additional options option: -o" $LINENO

#--------------------------------------------------------------------------------------------------------------------------------------------------

## Extra options
# other options
MERGE_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${OPTIONS}`
PRIORITIZATION_STRING_PARSED=`sed -e "s/'//g" <<< ${PRIORITIZATION_STRING}`


# Create output file name
OUTVCF=somaticvariants.vcf.gz


#--------------------------------------------------------------------------------------------------------------------------------------------------
## Merge VCFs with GATK's CombineVariants
#--------------------------------------------------------------------------------------------------------------------------------------------------
## Record start time
logInfo "[CombineVariants] START."

## Run CombineVariants
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. CombineVariants error. Check tool log ${TOOL_LOG}. " ' INT TERM EXIT
${JAVA}/java -jar ${GATK}/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R ${REFGEN} \
	--variant:strelka ${STRELKA_VCF} \
	--variant:mutect ${MUTECT_VCF} \
        -nt ${THR} \
        -genotypeMergeOptions PRIORITIZE -priority ${PRIORITIZATION_STRING_PARSED} \
	-o ${OUTVCF} \
	${MERGE_OPTIONS_PARSED}
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO

checkFile ${OUTVCF} "Failed to create output merged VCF file." $LINENO

#----------------------------------------------------------------------------------------------------------------------------------------------
## Post-Processing
#----------------------------------------------------------------------------------------------------------------------------------------------

## Check for the creation of the output VCF file
checkFile ${OUTVCF} "Output VCF is empty." $LINENO

## Open read permissions to the user group
chmod g+r ${OUTVCF}

## Check for the creation of the output VCF index file
checkFile ${OUTVCF}.tbi "Output VCF index is empty." $LINENO

## Open read permissions to the user group
chmod g+r ${OUTVCF}.tbi




#----------------------------------------------------------------------------------------------------------------------------------------------

logInfo "Merge VCFs completed. Find merged result in ${OUTVCF}."

#----------------------------------------------------------------------------------------------------------------------------------------------
