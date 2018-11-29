#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------
## strelka.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Script to run Strelka on a normal and tumor bam file
# 
#############################################################################
 USAGE:
 strelka.sh        
		   -s		<sample_name>
		   -B           <normal_bam> 
                   -T		<tumor_bam>
		   -g		<reference_genome_fasta>
                   -v           <outputVCF>
                   -I           <install_path>
		   -c		<python_config_file>
		   -t		<threads>
		   -F		<shared_functions>
	           -o		<additonal options>
		   -d		Turn on debug mode
                   -h           Display this usage/help text(No arg)
                   
EXAMPLES:
strelka.sh -h
strelka.sh -B normal -T tumor -g reference_genome -v outputVCF -I path/to/install -F path/to/shared_functions -o option

NOTES: 

#############################################################################


DOCS

set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=strelka.sh
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
while getopts ":hs:B:T:g:v:I:M:c:t:F:o:d" OPT
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
		I )  # Install path
			INSTALL=${OPTARG}
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
ERRLOG=${SAMPLE}.strelka.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.strelka.log

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

checkVar "${INSTALL+x}" "Missing install directory option: -I" $LINENO
checkDir ${INSTALL} "Reason= directory ${INSTALL} is not a directory or does not exist." $LINENO

checkVar "${BCF+x}" "Missing BCFTools directory option: -M" $LINENO
checkDir ${BCF} "Reason= BCFTools directory ${BCF} is not a directory or does not exist." $LINENO

checkVar "${CONFIG+x}" "Missing python configuration file: -c" $LINENO
#checkFile ${CONFIG} "Python configuration file ${CONFIG} is empty or does not exist." $LINENO

checkVar "${THREADS+x}" "Missing number of threads option: -t" $LINENO

checkVar "${SHARED_FUNCTIONS+x}" "Missing shared functions option: -F" $LINENO
checkFile ${SHARED_FUNCTIONS} "Shared functions file ${SHARED_FUNCTIONS} is empty or does not exist." $LINENO

checkVar "${OPTIONS}" "Missing additional options option: -o" $LINENO




#CHECKV_PATH=="`dirname "$0"`" #parse the directory of this function to locate the checkfiles script
#source ${CHECKV_PATH}/check_variable.sh

#CHECKF_PATH=="`dirname "$0"`" #parse the directory of this function to locate the checkfiles script
#source ${CHECKF_PATH}/check_file.sh



#--------------------------------------------------------------------------------------------------------------------------------------------------
STRELKA_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${OPTIONS}`





#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform Strelka variant calling
#--------------------------------------------------------------------------------------------------------------------------------------------------
## Record start time
logInfo "[Strelka] START."

## first configure the strelka run
#
## Add command trap here
#
${INSTALL}/bin/configureStrelkaSomaticWorkflow.py \
    --tumorBam=${TUMOR} \
    --normalBam=${NORMAL} \
    --referenceFasta=${REFGEN} \
    --runDir=./strelka \
    $STRELKA_OPTIONS_PARSED
#
## Add exitCode check
#

## Run the strelka workflow. We choose local so it will not spawn new jobs on the cluster
#
## Add command trap here
#
./strelka/runWorkflow.py -m local -j ${THR}

#
## Add exitCode check
#

checkFile ./strelka/results/variants/somatic.indels.vcf.gz "Output somatic indels file failed to create." $LINENO
checkFile ./strelka/results/variants/somatic.snvs.vcf.gz "Output somatic SNV file failed to create." $LINENO

#----------------------------------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------------------------------------
## Post-Processing
#----------------------------------------------------------------------------------------------------------------------------------------------

## Clean up strelka indel output
# zcat ./strelka/results/variants/somatic.indels.vcf.gz | ### Find fixStrelka_GT_indels.pl | tabix bgzip -f > somatic.indels.fixed.vcf.gz
#
## Check exitCode
#

## Clean up strelka snv output
# zcat ./strelka/results/variants/somatic.snvs.vcf.gz | ### Find fixStrelka_GT_snvs.pl | tabix bgzip -f > somatic.snvs.fixed.vcf.gz
#
## Check exitCode 
#

## Combine vcfs into single output
${BCF}/bcftools merge -m any -f PASS,. --force-sample ./strelka/results/variants/*.vcf.gz > ${SAMPLE}.vcf
gzip ${SAMPLE}.vcf.gz
#Not sure what this does -> | ${BCF}/bcftools plugin fill-AN-AC | ${BCF}/bcftools filter -i 'SUM(AC)>1' > ${SAMPLE}.vcf.gz

#----------------------------------------------------------------------------------------------------------------------------------------------

logInfo "Strelka workflow completed."

#----------------------------------------------------------------------------------------------------------------------------------------------
