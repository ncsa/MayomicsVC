#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## aggregate_gvcfs.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Aggregate gvcf files produced from the haplotypecaller for joint calling
# with GenotypeGVCFs
# 
#############################################################################

 USAGE:
aggregate_gvcfs.sh    -b           <sample1.g.vcf[,sample2.g.vcf,...]>
                      -S           </path/to/gatk/executable>
                      -e           </path/to/java_options_file>
                      -F           </path/to/shared_functions.sh>
                      -I           <genomic_intervals>
                      -d           turn on debug mode

 EXAMPLES:
aggregate_gvcfs.sh -h
aggregate_gvcfs.sh -b sample1.g.vcf,sample2.g.vcf,sample3.g.vcf -S /path/to/gatk/executable  -e /path/to/java_options_file -F /path/to/shared_functions.sh -I chr20 -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=aggregate_gvcfs.sh
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
while getopts ":hb:S:e:F:I:d" OPT
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
ERRLOG=gvcfs.aggregation.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=gvcfs.aggregate_gatk.log
truncate -s 0 ${TOOL_LOG}

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTGVCF+x}" "Missing input gvcf option: -b" $LINENO
for GVCF in $(echo ${INPUTGVCF} | sed "s/,/ /g")
do
        checkFile ${GVCF} "Input variants file ${GVCF} is empty or does not exist." $LINENO
        checkFile ${GVCF}.idx "Input variants index file ${GVCF}.idx is empty or does not exist." $LINENO
done

checkVar "${GATKEXE+x}" "Missing GATKEXE path option: -S" $LINENO
checkFileExe ${GATKEXE} "REASON=GATK file ${GATKEXE} is not an executable or does not exist." $LINENO

checkVar "${JAVA_OPTS_FILE+x}" "Missing file of JAVA path and options: -e" $LINENO
source ${JAVA_OPTS_FILE}
checkVar "${JAVA_PATH+x}" "Missing JAVA path from: -e ${JAVA_OPTS_FILE}" $LINENO
checkDir ${JAVA_PATH} "REASON=JAVA path ${JAVA_PATH} is not a directory or does not exist." $LINENO
checkVar "${JAVA_OPTS+x}" "Missing JAVA options from: -e ${JAVA_OPTS_FILE}" $LINENO

checkVar "${INTERVALS+x}" "Missing Intervals option: -I ${INTERVALS}" $LINENO




#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------






## Defining file names
VARIANTS=$( echo ${INPUTGVCF} | sed "s/,/ --variant /g" | tr "\n" " " )
OUTDB=variantsDB

rm -rf ${OUTDB}



#-------------------------------------------------------------------------------------------------------------------------------
## Perform Joint Calling with GATK GenotypeGVCFs
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Aggregating the per sample variants files"

## gatk genotypegvcf command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GenomicsDBImport aggregation error. " ' INT TERM EXIT
# Memory setting based on GATK's recommendation here: https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/JointGenotypingWf.wdl#L313
${GATKEXE} --java-options "-Xmx4g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path ${OUTDB} --intervals ${INTERVALS} --batch-size 50 --consolidate true --variant ${VARIANTS} --reader-threads 5 -ip 500 >> ${TOOL_LOG} 2>&1 
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] Aggregation of input GVCFs complete."




logInfo "[GATKEXE] Tarring the GVCFs Databse"
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Database tar error. " ' INT TERM EXIT
tar -cf ${OUTDB}.tar ${OUTDB}
EXITCODE=$?
trap - INT TERM EXIT
checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] GVCFs database tarring complete."

#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output Database and its index. Open read permissions to the user group
checkFile ${OUTDB}.tar "Output variants file ${OUTDB}.tar is empty." $LINENO

chmod g+r ${OUTDB}.tar


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
