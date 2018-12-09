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
                   -M           <BCFTools_path>
                   -I           <strelka_install_path>
                   -S           <Samtools_path>
                   -Z           <bgzip_path>
		   -t		<threads>
		   -F		<shared_functions>
	           -o		<additonal options>
		   -d		Turn on debug mode
                   -h           Display this usage/help text(No arg)
                   
EXAMPLES:
strelka.sh -h
strelka.sh -B normal.fastq -T tumor.fastq -g reference_genome.fasta -I /path/to/strelka/install -M /path/to/BCFTools -S /path/to/samtools -Z /path/to/bgzip -F /path/to/MayomicsVC/shared_functions.sh -o "'--extra_option'"

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
while getopts ":hs:B:T:g:I:M:S:Z:t:F:o:d" OPT
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
		I )  # Install path
			INSTALL=${OPTARG}
			checkArg
			;;
		M )  # BCF Tools path
			BCF=${OPTARG}
			checkArg
			;;
		S )  # Samtools path
			SAMTOOLS=${OPTARG}
			checkArg
			;;
		Z )  # bgzip path
			BGZIP=${OPTARG}
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
TOOL_LOG=${SAMPLE}.strelka_tool.log
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

#checkVar "${OUTVCF+x}" "Missing output VCF option: -v" $LINENO
#checkFile ${OUTVCF} "Output VCF file ${OUTVCF} is empty or does not exist." $LINENO

checkVar "${INSTALL+x}" "Missing install directory option: -I" $LINENO
checkDir ${INSTALL} "Reason= directory ${INSTALL} is not a directory or does not exist." $LINENO

checkVar "${BCF+x}" "Missing BCFTools directory option: -M" $LINENO
checkDir ${BCF} "Reason= BCFTools directory ${BCF} is not a directory or does not exist." $LINENO

checkVar "${SAMTOOLS+x}" "Missing Samtools directory option: -S" $LINENO
checkDir ${SAMTOOLS} "Reason= Samtools directory ${SAMTOOLS} is not a directory or does not exist." $LINENO

checkVar "${BGZIP+x}" "Missing bgzip directory option: -Z" $LINENO
checkDir ${BGZIP} "Reason= Bgzip directory ${BGZIP} is not a directory or does not exist." $LINENO

#checkVar "${CONFIG+x}" "Missing python configuration file: -c" $LINENO
#checkFile ${CONFIG} "Python configuration file ${CONFIG} is empty or does not exist." $LINENO

checkVar "${THR+x}" "Missing number of threads option: -t" $LINENO

checkVar "${SHARED_FUNCTIONS+x}" "Missing shared functions option: -F" $LINENO
checkFile ${SHARED_FUNCTIONS} "Shared functions file ${SHARED_FUNCTIONS} is empty or does not exist." $LINENO

checkVar "${OPTIONS}" "Missing additional options option: -o" $LINENO

#--------------------------------------------------------------------------------------------------------------------------------------------------

## Parsing extra options string
STRELKA_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${OPTIONS}`





#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform Strelka variant calling
#--------------------------------------------------------------------------------------------------------------------------------------------------
## Record start time
logInfo "[Strelka] START."

## first configure the strelka run
TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in configuring Strelka somatic workflow.  " ' INT TERM EXIT
${INSTALL}/configureStrelkaSomaticWorkflow.py \
    --tumorBam=${TUMOR} \
    --normalBam=${NORMAL} \
    --referenceFasta=${REFGEN} \
    --runDir=./strelka \
    $STRELKA_OPTIONS_PARSED >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO


## Run the strelka workflow. We choose local so it will not spawn new jobs on the cluster
TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in execution of runWorkflow.py. " ' INT TERM EXIT
./strelka/runWorkflow.py -m local -j ${THR} >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO

logInfo "[strelka] Finished running strelka workflow successfully for ${SAMPLE}"

checkFile ./strelka/results/variants/somatic.indels.vcf.gz "Output somatic indels file failed to create." $LINENO
checkFile ./strelka/results/variants/somatic.snvs.vcf.gz "Output somatic SNV file failed to create." $LINENO

#----------------------------------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------------------------------------
## Reformat SNV and Indel Genotypes
#----------------------------------------------------------------------------------------------------------------------------------------------

## Clean up strelka indel output
#zcat ./strelka/results/variants/somatic.indels.vcf.gz | perl ${FIX_INDEL_GT} | gzip somatic.indels.fixed.vcf.gz

#
## Check exitCode
#
EXITCODE=$?
#trap
checkExitcode ${EXITCODE} $LINENO

## Clean up strelka snv output
#zcat ./strelka/results/variants/somatic.snvs.vcf.gz | perl ${FIX_SNV_GT} | gzip somatic.snvs.fixed.vcf.gz

#
## Check exitCode 
#
EXITCODE=$?
#trap
checkExitcode ${EXITCODE} $LINENO

#----------------------------------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------------------------------------
## Post-Processing: VCF Merge
#----------------------------------------------------------------------------------------------------------------------------------------------

## Replace sample names
normal_sample_name=`${SAMTOOLS}/samtools view ${NORMAL} -H | awk '/@RG/ { for (i=1;i<=NF;i++) { if ($i ~ /SM:/) { sub("SM:","",$i); print $i; exit; } } }'`
tumor_sample_name=`${SAMTOOLS}/samtools view ${TUMOR} -H | awk '/@RG/ { for (i=1;i<=NF;i++) { if ($i ~ /SM:/) { sub("SM:","",$i); print $i; exit; } } }'`

## Combine vcfs into single output
logInfo "[BCFTools] Merging strelka output VCFs."

TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. BCFtools merge error. " ' INT TERM EXIT
${BCF}/bcftools merge -m any -f PASS,. --force-sample ./strelka/results/variants/*.vcf.gz > ${SAMPLE}.vcf 2>>${TOOL_LOG}
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

logInfo "[BCFTools] Finished VCF merge."

#----------------------------------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------------------------------------
## Post-Processing: bgzip and tabix
#----------------------------------------------------------------------------------------------------------------------------------------------

cat ${SAMPLE}.vcf | sed -e "/^#CHROM/ s/NORMAL/$normal_sample_name/g" -e "/^#CHROM/ s/TUMOR/$tumor_sample_name/g" > ${SAMPLE}.vcf.tmp 2>>${TOOL_LOG}




## BGZip merged vcf
logInfo "[bgzip] Zipping merged VCF."

TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. BGZIP error. " ' INT TERM EXIT
${BGZIP}/bgzip -c ${SAMPLE}.vcf.tmp > ${SAMPLE}.vcf.bgz 2>>${TOOL_LOG}
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

logInfo "[bgzip] Finished VCF zipping."




## Generate Tabix index
logInfo "[tabix] Creating Tabix index."

TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. BCFtools tabix error. " ' INT TERM EXIT
${BCF}/bcftools tabix -f -p vcf ${SAMPLE}.vcf.bgz >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

logInfo "[tabix] Finished index generation."

#Not sure what this does -> | ${BCF}/bcftools plugin fill-AN-AC | ${BCF}/bcftools filter -i 'SUM(AC)>1' > ${SAMPLE}.vcf.gz

#----------------------------------------------------------------------------------------------------------------------------------------------

logInfo "Strelka workflow completed."

#----------------------------------------------------------------------------------------------------------------------------------------------
