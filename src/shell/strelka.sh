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
		   -N           <normal_bam> 
                   -T		<tumor_bam>
		   -g		<reference_genome_fasta>
                   -B           <BCFTools_path>
                   -I           <strelka_install_path>
                   -S           <Samtools_path>
                   -Z           <bgzip_path>
		   -t		<threads>
                   -e           <environmental_profile>
		   -F		<shared_functions>
	           -i           <indel_GT_fix_perl_script>
                   -p           <snp_GT_fix_perl_script>
                   -o		<additonal strelka options>
                   -O		<additonal bcf options>
		   -d		Turn on debug mode
                   -h           Display this usage/help text(No arg)
                   
EXAMPLES:
strelka.sh -h
strelka.sh -N normal.fastq -T tumor.fastq -g reference_genome.fasta -I /path/to/strelka/install -B /path/to/BCFTools -S /path/to/samtools -Z /path/to/bgzip -e /path/to/envprofile.file -F /path/to/MayomicsVC/shared_functions.sh -i /path/to/fix_indels.pl -p /path/to/fix_snps.pl -o "'--extra_option'" -O "'extra_bcf_options'"

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
while getopts ":hs:N:T:g:I:B:S:Z:t:e:F:i:p:o:d" OPT
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
                N )  # Normal sample BAM
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
		B )  # BCF Tools path
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
                e )  # Path to file with environmental profile variables
                        ENV_PROFILE=${OPTARG}
                        checkArg
                        ;;
		F )  # Shared functions 
                        SHARED_FUNCTIONS=${OPTARG}
                        checkArg
                        ;;
		i )  # Path to fixStrelka_GT_indels.pl
			FIX_INDEL_GT=${OPTARG}
			checkArg
			;;
		p )  # Path to fixStrelka_GT_snvs.pl
			FIX_SNV_GT=${OPTARG}
			checkArg
			;;
		o )  # Extra strelka options
                        STRELKA_OPTIONS=${OPTARG}
                        checkArg
                        ;;
                O )  # Extra bcf tools options
                        BCFTOOLS_OPTIONS=${OPTARG}
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

source "${SHARED_FUNCTIONS}"

## Check if Sample Name variable exists
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Create log for JOB_ID/script and tool
ERRLOG=${SAMPLE}.strelka_variant_calling.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=${SAMPLE}.strelka.log                                                                                                            
truncate -s 0 ${TOOL_LOG}

## Send manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with environmental profile variables
checkVar "${ENV_PROFILE+x}" "Missing environmental profile option: -e" $LINENO
source ${ENV_PROFILE}


## Check if input files, directories, and variables are non-zero
checkVar "${NORMAL+x}" "Missing normal BAM option: -B" $LINENO
checkFile ${NORMAL} "Input normal BAM file ${NORMAL} is empty or does not exist." $LINENO

checkVar "${TUMOR+x}" "Missing tumor BAM option: -T" $LINENO
checkFile ${TUMOR} "Input tumor BAM file ${TUMOR} is empty or does not exist." $LINENO

checkVar "${REFGEN+x}" "Missing reference genome option: -g" $LINENO
checkFile ${REFGEN} "Input tumor BAM file ${REFGEN} is empty or does not exist." $LINENO

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

checkVar "${FIX_INDEL_GT+x}" "Missing Fix GT Indel perl script option: -i" $LINENO
checkFile ${FIX_INDEL_GT} "Fix GT INDEL perl script ${FIX_INDEL_GT} is empty or does not exist." $LINENO

checkVar "${FIX_SNV_GT+x}" "Missing Fix GT SNV perl script option: -p" $LINENO
checkFile ${FIX_SNV_GT} "Fix GT SNV perl script ${FIX_SNV_GT} is empty or does not exist." $LINENO

checkVar "${STRELKA_OPTIONS}" "Missing additional strelka options option: -o" $LINENO
checkVar "${BCFTOOLS_OPTIONS}" "Missing additional bcftools options option: -o" $LINENO

#--------------------------------------------------------------------------------------------------------------------------------------------------

## Parsing extra options string
STRELKA_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${STRELKA_OPTIONS}`
BCFTOOLS_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${BCFTOOLS_OPTIONS}`





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
TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in execution of fix Indel GT script. " ' INT TERM EXIT
zcat ./strelka/results/variants/somatic.indels.vcf.gz | perl ${FIX_INDEL_GT} > somatic.indels.fixed.vcf 2>> ${TOOL_LOG} 
${BGZIP}/bgzip -c somatic.indels.fixed.vcf > somatic.indels.fixed.vcf.bgz 2>> ${TOOL_LOG}
${BCF}/bcftools tabix -f -p vcf somatic.indels.fixed.vcf.bgz >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO

## Clean up strelka snv output
TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in execution fix SNV GT script. " ' INT TERM EXIT
zcat ./strelka/results/variants/somatic.snvs.vcf.gz | perl ${FIX_SNV_GT} > somatic.snvs.fixed.vcf 2>> ${TOOL_LOG} 
${BGZIP}/bgzip -c somatic.snvs.fixed.vcf > somatic.snvs.fixed.vcf.bgz 2>> ${TOOL_LOG}
${BCF}/bcftools tabix -f -p vcf somatic.snvs.fixed.vcf.bgz >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT

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
${BCF}/bcftools merge ${BCFTOOLS_OPTIONS_PARSED} ./*fixed.vcf.bgz > ${SAMPLE}.vcf 2>> ${TOOL_LOG}
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO

logInfo "[BCFTools] Finished VCF merge."

#----------------------------------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------------------------------------
## Post-Processing: bgzip and tabix
#----------------------------------------------------------------------------------------------------------------------------------------------

cat ${SAMPLE}.vcf | sed -e "/^#CHROM/ s/NORMAL/$normal_sample_name/g" -e "/^#CHROM/ s/TUMOR/$tumor_sample_name/g" > ${SAMPLE}.vcf.tmp 2>> ${TOOL_LOG}




## BGZip merged vcf
logInfo "[bgzip] Zipping merged VCF."

TRAP_LINE=$(($LINENO+1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. BGZIP error. " ' INT TERM EXIT
${BGZIP}/bgzip -c ${SAMPLE}.vcf.tmp > ${SAMPLE}.vcf.bgz 2>> ${TOOL_LOG}
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
