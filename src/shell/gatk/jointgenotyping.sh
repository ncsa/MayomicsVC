#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## jointgenotyping.sh MANIFEST, USAGE DOCS, SET CHECKS
#-------------------------------------------------------------------------------------------------------------------------------

read -r -d '' MANIFEST << MANIFEST

*****************************************************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*****************************************************************************

MANIFEST
echo -e "${MANIFEST}"








read -r -d '' DOCS << DOCS

#############################################################################
#
# Jointly call from a gvcf produced from the haplotypecaller/CombineGVCFs or
# a GenomicsDB workspace
# 
#############################################################################

 USAGE:
 jointgenotyping.sh   -b           <input.g.vcf or genomicsDB.tar>
                      -S           </path/to/gatk/executable>
                      -G           <reference_genome>
                      -D           <dbsnp.vcf>
                      -o           <extra_genotypegvcf_options>
                      -e           </path/to/java_options_file>
                      -F           </path/to/shared_functions.sh>
                      -I           <genomic_intervals>
                      -d           turn on debug mode

 EXAMPLES:
 jointgenotyping.sh -h
 jointgenotyping.sh -b input.g.vcf -S /path/to/gatk/executable -G reference.fa -D dbsnp.vcf -o "'--sample_ploidy 2 --useNewAFCalculator'"  -e /path/to/java_options_file -F /path/to/shared_functions.sh -I chr20 -d

 NOTE: In order for getops to read in a string arguments for -o (extra_genotypegvcf_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.
#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=jointgenotyping.sh
SGE_JOB_ID=TBD  # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID




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
## GETOPTS ARGUMENT PARSER
#-------------------------------------------------------------------------------------------------------------------------------

## Check if no arguments were passed
if (($# == 0))
then
        echo -e "\nNo arguments passed.\n\n${DOCS}\n"
        exit 1
fi

## Input and Output parameters
while getopts ":hb:S:G:D:o:e:F:I:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo -e "\n${DOCS}\n"
                        exit 0
                        ;;
                b )  # Full path to the input gvcfs or list of gvcfs
                        INPUTGVCF=${OPTARG}
                        checkArg
                        ;;
                S )  # Full path to gatk executable
                        GATKEXE=${OPTARG}
                        checkArg
                        ;;
                G )  # Full path to referance genome fasta file
                        REF=${OPTARG} 
                        checkArg
                        ;;
                D ) # Full path to DBSNP file
                        DBSNP=${OPTARG}
                        checkArg
                        ;;
                o ) # Extra options and arguments to haplotyper, input as a long string, can be empty if desired
                        GENOTYPEGVCF_OPTIONS=${OPTARG}
                        checkArg
                        ;;
                e )  # Path to file containing JAVA options to pass into the gatk command 
                        JAVA_OPTS_FILE=${OPTARG}
                        checkArg
                        ;;
                F )  # Path to shared_functions.sh
                        SHARED_FUNCTIONS=${OPTARG}
                        checkArg
                        ;;
                I )  # Genomic intervals overwhich to operate
                        INTERVALS=${OPTARG}
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






#-------------------------------------------------------------------------------------------------------------------------------
## PRECHECK FOR INPUTS AND OPTIONS
#-------------------------------------------------------------------------------------------------------------------------------


source ${SHARED_FUNCTIONS}

## Create log for JOB_ID/script
ERRLOG=${INTERVALS}.jointgenotyping.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=${INTERVALS}.jointgenotyping_gatk.log
truncate -s 0 ${TOOL_LOG}

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTGVCF+x}" "Missing input gvcf option: -b" $LINENO
checkFile ${INPUTGVCF} "Input variants file (or database) ${INPUTGVCF} is empty or does not exist." $LINENO
VARIANTS_OPTION="--variant ${INPUTGVCF}"
if [[ ! -z ${INPUTGVCF+x} ]] # if a tarred database file is input 
then
    tar -xf ${INPUTGVCF}
    GVCFDB=$(basename ${INPUTGVCF} .tar)
    checkDir ${GVCFDB} "Input variants database ${INPUTGVCF} is empty or malformated." $LINENO
    VARIANTS_OPTION="--variant gendb://${GVCFDB}"
fi

checkVar "${GATKEXE+x}" "Missing GATKEXE path option: -S" $LINENO
checkFileExe ${GATKEXE} "REASON=GATK file ${GATKEXE} is not an executable or does not exist." $LINENO

checkVar "${JAVA_OPTS_FILE+x}" "Missing file of JAVA path and options: -e" $LINENO
source ${JAVA_OPTS_FILE}
checkVar "${JAVA_PATH+x}" "Missing JAVA path from: -e ${JAVA_OPTS_FILE}" $LINENO
checkDir ${JAVA_PATH} "REASON=JAVA path ${JAVA_PATH} is not a directory or does not exist." $LINENO
checkVar "${JAVA_OPTS+x}" "Missing JAVA options from: -e ${JAVA_OPTS_FILE}" $LINENO

checkVar "${REF+x}" "Missing reference genome option: -G" $LINENO
checkFile ${REF} "Reference genome file ${REF} is empty or does not exist." $LINENO

checkVar "${DBSNP+x}" "Missing dbSNP option: -D" $LINENO
checkFile ${DBSNP} "DBSNP ${DBSNP} is empty or does not exist." $LINENO


checkVar "${GENOTYPEGVCF_OPTIONS+x}" "Missing extra haplotyper options option: -o" $LINENO

checkVar "${INTERVALS+x}" "Missing Intervals option: -I ${INTERVALS}" $LINENO



#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------
GENOTYPEGVCF_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${GENOTYPEGVCF_OPTIONS}`






## Defining file names
OUTVCF=${INTERVALS}.vcf






#-------------------------------------------------------------------------------------------------------------------------------
## Perform Joint Calling with GATK GenotypeGVCFs
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Joint calling via GenotypeGVCFs"

## gatk genotypegvcf command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GenotypeGVCFs joint calling error. " ' INT TERM EXIT
${GATKEXE} ${JAVA_OPTS} GenotypeGVCFs --reference ${REF} --output ${OUTVCF} ${VARIANTS_OPTION} --dbsnp ${DBSNP} ${GENOTYPEGVCF_OPTIONS_PARSED}--intervals ${INTERVALS}  >> ${TOOL_LOG} 2>&1 
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] Joint calling via GenotypeGVCFs complete."




#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output VCF and its index. Open read permissions to the user group
checkFile ${OUTVCF} "Output variants file ${OUTVCF} is empty." $LINENO
checkFile ${OUTVCF}.idx "Output variants index file is empty." $LINENO

chmod g+r ${OUTVCF}
chmod g+r ${OUTVCF}.idx


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
