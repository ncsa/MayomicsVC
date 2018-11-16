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
		   -S		<shared_functions_path>
	           -o		<additonal options>
		   -d		Turn on debug mode
                   -h           Display this usage/help text(No arg)
                   
EXAMPLES:
strelka.sh -h
strelka.sh -B normal -T tumor -g reference_genome -v outputVCF -I path/to/install -s path/to/shared_functions -o option

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


LOG_PATH="`dirname "$0"`"  ## Parse the directory of this script to locate the logging function script
source ${LOG_PATH}/log_functions.sh

#-------------------------------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------------------------------
## GETOPTS ARGUMENT PARSER
#-------------------------------------------------------------------------------------------------------------------------------

## Check if no arguments were passed
if (($# == 0))
then
	echo -e "\nNo arguments passed.\n\n${DOCS}\n"
	exit 1
fi

## Input and Output parameters
while getopts ":hs:s:B:T:g:v:I:d:o" OPT
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
	        c )  # Python configuration file
                        CONFIG=${OPTARG}
                        checkArg
                        ;;
	        t )  # Number of threads
                        THR=${OPTARG}
                        checkArg
                        ;;
		s )  # Shared functions 
                        FUNCTIONS=${OPTARG}
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
# Check if sample name is set
if [[ -z ${SAMPLE+x} ]] ## NOTE: ${VAR+x} is used for variable expansions, preventing unset variable error from set -o nounset. When $VAR is not set, we set it to "x" and throw the error.
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Missing tumor name option: -s"
        exit 1
fi


## Send Manifest to log
ERRLOG=${SAMPLE}.strelka.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${TUMOR}.strelka.log                 

echo "${MANIFEST}" >> "${ERRLOG}"
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if the sample name option was passed in.
if [[ -z ${SAMPLE} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing sample name option: -s"
fi

## Check if the sample name is present.
#if [[ ! -s ${SAMPLE} ]]
#then
#        EXITCODE=1
#        logError "$0 stopped at line $LINENO. \nREASON=SAMPLE NAME ${SAMPLE} is empty or does not exist."
#fi

## Check if the NOMRAL input BAM option was passed in.
if [[ -z ${NORMAL} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing NORMAL BAM option: -B"
fi

## Check if the NORMAL input BAM is present.
if [[ ! -s ${NORMAL} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=NORMAL BAM ${NORMAL} is empty or does not exist."
fi
## Check if the TUMOR input BAM option was passed in.
if [[ -z ${TUMOR} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing TUMOR BAM option: -T"
fi

## Check if the TUMOR input BAM is present.
if [[ ! -s ${TUMOR} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=TUMOR BAM ${TUMOR} is empty or does not exist."
fi

## Check if the reference file option was passed in.
if [[ -z ${REFGEN} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing reference genome option: -g"
fi

## Check if the reference file is present.
if [[ ! -s ${REFGEN} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome fasta file ${REFGEN} is empty or does not exist."
fi

## Check if the Output VCF option was passed in.
if [[ -z ${OUTVCF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing output VCF option: -o"
fi

## Check if the OUTPUT VCF file  is present.
if [[ ! -s ${OUTVCF} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome fasta file ${OUTVCF} is empty or does not exist."
fi

## Check if the Strelka INSTALL file option was passed in.
if [[ -z ${INSTALL} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Strelka INSTALL file option: -I"
fi

## Check if the Strelka INSTALL file is present.
if [[ ! -s ${INSTALL} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Strelka directory ${INSTALL} is not a directory or does not exist."
fi

## Check if the python configuration file option was passed in.
if [[ -z ${CONFIG} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing python configuration file option: -c"
fi

## Check if the python configuration file is present.
if [[ ! -s ${CONFIG} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Strelka directory ${CONFIG} is not a directory or does not exist."
fi

## Check if the number of threads option was passed in.
if [[ -z ${THR} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing number of threads option: -t"
fi

## Check if the number of threads is present.
if [[ ! -s ${THR} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Strelka directory ${THR} is not a directory or does not exist."
fi

## Check if the shared functions option was passed in.
if [[ -z ${FUNCTIONS} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Strelka INSTALL file option: -I"
fi

## Check if the shared functions file is present.
if [[ ! -s ${FUNCTIONS} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Strelka directory ${INSTALL} is not a directory or does not exist."
fi

## Check if the Extra Options option was passed in
if [[ -z ${OPTIONS+x} ]]
then
        EXITCODE=1
        logError "$0 stopped at line ${LINENO}. \nREASON=Missing Extra Options option: -o"
fi
## Check if the Extra Options is present.
if [[ ! -s ${OPTIONS} ]]
then
        EXITCODE=1
        logError "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REFGEN} is empty or does not exist."



#--------------------------------------------------------------------------------------------------------------------------------------------------
STRELKA_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${OPTIONS}`















#-----------------------------------------------------------------------------------
## debug mode
if [ "$DEBUG_MODE" == "YES" ];then set -x; fi
log "Command to recreate: $0 $@"

if [ ! -s $normal_bam ]
then
    $WORKFLOW_PATH/email.sh -f $normal_bam -m strelka.sh -p "$email_param" -l $LINENO
    exit 100;
else
    validate_bam $SAMTOOLS $normal_bam $output strelka
    if [ $? -gt 0 ]
    then
        $WORKFLOW_PATH/email.sh -f $normal_bam -m strelka.sh -p "$email_param" -l $LINENO
        exit 100;
    fi
fi

if [ ! -s $tumor_bam ]
then
    $WORKFLOW_PATH/email.sh -f $tumor_bam -m strelka.sh -p "$email_param" -l $LINENO
    exit 100;
else
    validate_bam $SAMTOOLS $tumor_bam $output strelka
    if [ $? -gt 0 ]
    then
        $WORKFLOW_PATH/email.sh -f $tumor_bam -m strelka.sh -p "$email_param" -l $LINENO
        exit 100;
    fi
fi

if [ ! -d $output/strelka ]
then
    mkdir -p $output/strelka
else
    rm -rv $output/strelka
    mkdir -p $output/strelka
fi
#----------------------------------------------------------------------------------------








#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform Strelka
#--------------------------------------------------------------------------------------------------------------------------------------------------
## Record start time
logInfo "[Strelka] START."

## first configure the strelka run
$STRELKA/bin/configureStrelkaSomaticWorkflow.py \
    --tumorBam=$TUMOR \
    --normalBam=$NORMAL \
    --referenceFasta=$REFGEN \
    --runDir=$output/strelka/ \
    $STRELKA_params
if [ $? -ne 0 ]
then
  $WORKFLOW_PATH/email.sh -f $output/$vcf -m somaticvariants.sh -s $STRELKA/bin/configureStrelkaSomaticWorkflow.py -p "email_param" -l $LINENO
  exit 100;   
fi

## then run the strelka workflow
## we choose local so it will not spawn new jobs on the cluster
$output/strelka/runWorkflow.py -m local -j $THR
if [ $? -ne 0 ]
then
  $WORKFLOW_PATH/email.sh -f $output/$vcf -m somaticvariants.sh -s $output/strelka/runWorkflow.py -p "email_param" -l $LINENO
  exit 100;   
fi

if [ ! -s  $output/strelka/results/variants/somatic.indels.vcf.gz ]
then
    $WORKFLOW_PATH/email.sh -f $output/strelka/results/variants/somatic.indels.vcf.gz -m somaticvariants.sh -s strelka.sh -p "$email_param" -l $LINENO -M "failed to create"
    exit 100;
fi
#----------------------------------------------------------------------------------------------------------------------------------------------











#----------------------------------------------------------------------------------------------------------------------------------------------
##POST-PROCESSING
#----------------------------------------------------------------------------------------------------------------------------------------------
# Clean up strelka output
zcat $output/strelka/results/variants/somatic.indels.vcf.gz | $PERL $WORKFLOW_PATH/fixStrelka_GT_indels.pl | $TABIX/bgzip -f > $output/strelka/results/variants/somatic.indels.fixed.vcf.gz

if [ ! -s  $output/strelka/results/variants/somatic.indels.fixed.vcf.gz ]
then
    $WORKFLOW_PATH/email.sh -f $output/strelka/results/variants/somatic.indels.fixed.vcf.gz -m somaticvariants.sh -s strelka.sh -p "$email_param" -l $LINENO -M "failed to create"
    exit 100;
fi

# Since the output is valid, move the file to its final destination
mv $output/strelka/results/variants/somatic.indels.fixed.vcf.gz $output/strelka/results/variants/somatic.indels.vcf.gz
$TABIX/tabix -f -p vcf $output/strelka/results/variants/somatic.indels.vcf.gz

if [ ! -s  $output/strelka/results/variants/somatic.snvs.vcf.gz ]
then
    $WORKFLOW_PATH/email.sh -f $output/strelka/results/variants/somatic.snvs.vcf.gz -m somaticvariants.sh -s strelka.sh -p "$email_param" -l $LINENO -M "failed to create"
    exit 1;
fi

# Clean up strelka output
zcat $output/strelka/results/variants/somatic.snvs.vcf.gz | $PERL $WORKFLOW_PATH/fixStrelka_GT_snvs.pl | $TABIX/bgzip -f > $output/strelka/results/variants/somatic.snvs.fixed.vcf.gz

if [ ! -s  $output/strelka/results/variants/somatic.snvs.fixed.vcf.gz ]
then
    $WORKFLOW_PATH/email.sh -f $output/strelka/results/variants/somatic.snvs.fixed.vcf.gz -m somaticvariants.sh -s strelka.sh -p "$email_param" -l $LINENO -M "failed to create"
    exit 100;
fi

# Since the output is valid, move the file to its final destination
mv $output/strelka/results/variants/somatic.snvs.fixed.vcf.gz $output/strelka/results/variants/somatic.snvs.vcf.gz
$TABIX/tabix -f -p vcf $output/strelka/results/variants/somatic.snvs.vcf.gz

## combine vcfs into a single output
if [ -f $output/$vcf ]; then rm $output/$vcf; fi

$WORKFLOW_PATH/combinevcf.sh \
    -i $output/strelka/results/variants/somatic.indels.vcf.gz \
    -i $output/strelka/results/variants/somatic.snvs.vcf.gz \
    -o $output/$vcf \
    -t $TOOL_INFO \
    -m $MEMORY_INFO \
    -b $tumor_bam \
    -a
if [ $? -ne 0 ]
then
  $WORKFLOW_PATH/email.sh -f $output/$vcf -m somaticvariants.sh -s combinevcf.sh -p "email_param" -l $LINENO
  exit 100;   
fi

normal_sample_name=`$SAMTOOLS/samtools view $normal_bam -H | awk '/@RG/ { for (i=1;i<=NF;i++) { if ($i ~ /SM:/) { sub("SM:","",$i); print $i; exit; } } }'`
tumor_sample_name=`$SAMTOOLS/samtools view $tumor_bam -H | awk '/@RG/ { for (i=1;i<=NF;i++) { if ($i ~ /SM:/) { sub("SM:","",$i); print $i; exit; } } }'`

# replace NORMAL and TUMOR in the VCF column header line with the actual sample names
if [ -f $output/$vcf.tmp ]; then rm $output/$vcf.tmp; fi

zcat $output/$vcf | sed -e "/^#CHROM/ s/NORMAL/$normal_sample_name/g" -e "/^#CHROM/ s/TUMOR/$tumor_sample_name/g" > $output/$vcf.tmp

$TABIX/bgzip -c $output/$vcf.tmp > $output/$vcf
rm -v $output/$vcf.tmp

$TABIX/tabix -f -p vcf $output/$vcf

END=$(date +%s)
DIFF=$(( $END - $START ))
log "strelka.sh to create $output/$vcf took $DIFF seconds"
