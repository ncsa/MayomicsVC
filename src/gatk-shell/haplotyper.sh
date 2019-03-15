#!/usr/bin/env bash

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


###########################################################################################################################
#
# Perform GATK's HaplotypeCaller variant caller on the bam produced in the Deduplication stage of the Mayomics workflow.
# bqsr.sh must be run before to this stage to calculate the required modification of the quality scores.
# Step 2/3 in Single Sample Variant Calling.
#
###########################################################################################################################

 USAGE:
 haplotyper.sh     -s 	<sample_name>
                   -S	</path/to/gatk/executable>
                   -G	<reference_genome>
                   -t	<threads>
                   -b	<sorted.deduped.bam>
                   -D	<dbsnp.vcf>
                   -o	<extra_haplotyper_options>
                   -e   </path/to/java_options_file>
                   -F   </path/to/shared_functions.sh>
                   -V   VCF Source Field (default: Sentieon)
                   -d   turn on debug mode

 EXAMPLES:
 haplotyper.sh -h
 haplotyper.sh -s sample -S /path/to/gatk/executable -G reference.fa -t 12 -b sorted.deduped.bam -D dbsnp.vcf -o "'--sample-ploidy 2 -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest '" -e /path/to/java_options_file -F </path/to/shared_functions.sh> -d 

NOTE: In order for getops to read in a string arguments for -o (extra_haplotyper_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.
###########################################################################################################################


DOCS

set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=haplotyper.sh
SGE_JOB_ID=TBD  # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID
SOURCE_FIELD="Sentieon"


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

while getopts ":hs:S:G:t:b:D:o:e:F:dV:" OPT
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
		o ) #Extra options and arguments to haplotyper, input as a long string, can be empty if desired
			HAPLOTYPER_OPTIONS=${OPTARG}
			checkArg
			;;
        e )  # Path to file containing JAVA options to pass into the gatk command 
            JAVA_OPTS_FILE=${OPTARG}
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
	    V ) # Set the source field overriding 'Sentieon'
	        SOURCE_FIELD=${OPTARG}
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
truncate -s 0 ${SAMPLE}.haplotype_gatk.log

echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with java path and options variables
checkVar "${JAVA_OPTS_FILE+x}" "Missing file of JAVA path and options: -e" $LINENO
source ${JAVA_OPTS_FILE}
checkVar "${JAVA_PATH+x}" "Missing JAVA path from: -e ${JAVA_OPTS_FILE}" $LINENO
checkDir ${JAVA_PATH} "REASON=JAVA path ${JAVA_PATH} is not a directory or does not exist." $LINENO
checkVar "${JAVA_OPTS+x}" "Missing JAVA options from: -e ${JAVA_OPTS_FILE}" $LINENO

checkVar "${GATKEXE+x}" "Missing GATK path option: -S" $LINENO
checkFileExe ${GATKEXE} "REASON=GATK file ${GATKEXE} is not executable or does not exist." $LINENO

checkVar "${NTHREADS+x}" "Missing threads option: -t" $LINENO

checkVar "${REF+x}" "Missing reference genome option: -G" $LINENO
checkFile ${REF} "Reference genome file ${REF} is empty or does not exist." $LINENO

checkVar "${INPUTBAM+x}" "Missing input BAM option: -b" $LINENO
checkFile ${INPUTBAM} "Input BAM ${INPUTBAM} is empty or does not exist." $LINENO
checkFile ${INPUTBAM}.bai "Input BAM index ${INPUTBAM} is empty or does not exist." $LINENO

checkVar "${DBSNP+x}" "Missing dbSNP option: -D" $LINENO
checkFile ${DBSNP} "DBSNP ${DBSNP} is empty or does not exist." $LINENO


checkVar "${HAPLOTYPER_OPTIONS+x}" "Missing extra haplotyper options option: -o" $LINENO


#--------------------------------------------------------------------------------------------------------------------------
HAPLOTYPER_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${HAPLOTYPER_OPTIONS}`

TOOL_LOG=${SAMPLE}.haplotype_gatk.log








#-------------------------------------------------------------------------------------------------------------------------
## Perform HaplotypeCaller with GATK.
#-------------------------------------------------------------------------------------------------------------------------


## Record start time
logInfo "[HaplotypeCaller] START."


#Execute GATK with the HaplotypeCaller algorithm
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Error in GATK HaplotypeCaller. " ' INT TERM EXIT
${GATKEXE} HaplotypeCaller --native-pair-hmm-threads ${NTHREADS} --reference ${REF} --input ${INPUTBAM} --output ${SAMPLE}.g.vcf --dbsnp ${DBSNP} ${HAPLOTYPER_OPTIONS_PARSED} --emit-ref-confidence GVCF >> ${TOOL_LOG} 2>&1

EXITCODE=$?
trap - INT TERM EXIT


checkExitcode ${EXITCODE} $LINENO
logInfo "[HaplotypeCaller] Finished running successfully. Output: ${SAMPLE}.g.vcf"
#-------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------


## Check for the creation of the output VCF file
checkFile ${SAMPLE}.g.vcf "Output VCF is empty." $LINENO

## Open read permissions to the user group
chmod g+r ${SAMPLE}.g.vcf

## Check for the creation of the output VCF index file 
checkFile ${SAMPLE}.g.vcf.idx "Output VCF index is empty." $LINENO

## Open read permissions to the user group
chmod g+r ${SAMPLE}.g.vcf.idx


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
