#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------------
## haplotyper.sh MANIFEST, USAGE DOCS, SET CHECKS
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


#######################################################################################################################################################
#
# Perform Sentieon's Haplotyper variant caller on the bam produced in the Deduplication stage of the Mayomics workflow.
# bqsr.sh must be run before to this stage to calculate the required modification of the quality scores.
# Step 2/3 in Single Sample Variant Calling.
#
########################################################################################################################################################

 USAGE:
 Haplotyper.sh     -s 	<sample_name>
		   -S	</path/to/sentieon>
		   -G	<reference_genome>
		   -t	<threads>
		   -b	<sorted.deduped.realigned.bam>
		   -D	<dbsnp.vcf>
		   -r	<recal_data.table>
		   -o	<extra_haplotyper_options>
                   -e   </path/to/env_profile_file>
                   -F   </path/to/shared_functions.sh>
		   -d   turn on debug mode

 EXAMPLES:
 Haplotyper.sh -h
 Haplotyper.sh -s sample -S /path/to/sentieon_directory -G reference.fa -t 12 -b sorted.deduped.realigned.recalibrated.bam -D dbsnp.vcf -r recal_data.table -o "'--emit_mode variant --gq_bands 1-60,60-99/19,99 --min_base_qual 10 --pcr_indel_model CONSERVATIVE --phasing 1 --ploidy 2 --prune_factor 2'" -e /path/to/env_profile_file -F </path/to/shared_functions.sh> -d 

NOTE: In order for getops to read in a string arguments for -o (extra_haplotyper_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.
##########################################################################################################################################################


DOCS

set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=Haplotyper.sh
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

while getopts ":hs:S:G:t:b:D:r:o:e:F:d" OPT
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
		S ) # Full path to sentieon directory
			SENTIEON=${OPTARG}
			checkArg
			;;
		G ) # Full path to referance genome fasta file
			REF=${OPTARG}
			checkArg
			;;
		t ) # Number of threads available
			NTHREADS=${OPTARG}
			checkArg
			;;
		b ) # Full path to BAM used as input
			INPUTBAM=${OPTARG}
			checkArg
			;;
		D ) # Full path to DBSNP file
			DBSNP=${OPTARG}
			checkArg
			;;
		r ) #Full path to the recal_data.table created in the BQSR step
			RECAL=${OPTARG}
			checkArg
			;;
		o ) #Extra options and arguments to haplotyper, input as a long string, can be empty if desired
			HAPLOTYPER_OPTIONS=${OPTARG}
			checkArg
			;;
                e )  # Path to file with environmental profile variables
                        ENV_PROFILE=${OPTARG}
                        checkArg
                        ;;
		F ) # Path to shared_functions.sh
			SHARED_FUNCTIONS=${OPTARG}
			checkArg
			;;
		d ) # Turn on debug mode. Turn on debug mode. Initiates 'set -x' to print all text. Invoked with -d
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


source ${SHARED_FUNCTIONS}


## Check if sample name is set
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Send Manifest to log
ERRLOG=${SAMPLE}.haplotyper.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.haplotype_sentieon.log

echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with environmental profile variables
checkVar "${ENV_PROFILE+x}" "Missing environmental profile option: -e" $LINENO
source ${ENV_PROFILE}

checkVar "${SENTIEON+x}" "Missing Sentieon path option: -S" $LINENO
checkDir ${SENTIEON} "REASON=Sentieon directory ${SENTIEON} is not a directory or does not exist." $LINENO

checkVar "${NTHREADS+x}" "Missing threads option: -t" $LINENO

checkVar "${REF+x}" "Missing reference genome option: -G" $LINENO
checkFile ${REF} "Reference genome file ${REF} is empty or does not exist." $LINENO

checkVar "${INPUTBAM+x}" "Missing input BAM option: -b" $LINENO
checkFile ${INPUTBAM} "Input BAM ${INPUTBAM} is empty or does not exist." $LINENO
checkFile ${INPUTBAM}.bai "Input BAM index ${INPUTBAM} is empty or does not exist." $LINENO

checkVar "${DBSNP+x}" "Missing dbSNP option: -D" $LINENO
checkFile ${DBSNP} "DBSNP ${DBSNP} is empty or does not exist." $LINENO


RECAL_OPTION=""
if [[ ! -z ${RECAL+x} ]] # if unset return null, if null or set, return x
then
	## Check if the Recal_data.table file produced in BQSR is present
	if [[ ! -s ${RECAL} ]]
	then
		EXITCODE=1
		logError "$0 stopped at line $LINENO. \nREASON=RECAL_DATA.TABLE ${RECAL} is empty or does not exist." 
	fi
RECAL_OPTION="-q ${RECAL}"
fi

#checkVar "${RECAL+x}" "Missing RECAL_DATA.TABLE option: -r" $LINENO
#checkFile ${RECAL} "Recal table from BQSR ${RECAL} is empty or does not exist." $LINENO

checkVar "${HAPLOTYPER_OPTIONS+x}" "Missing extra haplotyper options option: -o" $LINENO


#--------------------------------------------------------------------------------------------------------------------------------------------------
HAPLOTYPER_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${HAPLOTYPER_OPTIONS}`










#--------------------------------------------------------------------------------------------------------------------------------------------------
## Perform Haplotyper with Sentieon.
#--------------------------------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[Haplotyper] START."


#Execute Sentieon with the Haplotyper algorithm
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in Sentieon Haplotyper. " ' INT TERM EXIT
${SENTIEON}/bin/sentieon driver -t ${NTHREADS} -r ${REF} -i ${INPUTBAM} ${RECAL_OPTION} --algo Haplotyper ${HAPLOTYPER_OPTIONS_PARSED} -d ${DBSNP} ${SAMPLE}.vcf >> ${SAMPLE}.haplotype_sentieon.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT


checkExitcode ${EXITCODE} $LINENO
logInfo "[Haplotyper] Finished running successfully. Output: ${SAMPLE}.vcf"
#------------------------------------------------------------------------------------------------------------------------------------







#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for the creation of the output VCF file
checkFile ${SAMPLE}.vcf "Output VCF is empty." $LINENO

## Open read permissions to the user group
chmod g+r ${SAMPLE}.vcf

## Check for the creation of the output VCF index file 
checkFile ${SAMPLE}.vcf.idx "Output VCF index is empty." $LINENO

## Open read permissions to the user group
chmod g+r ${SAMPLE}.vcf.idx


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
