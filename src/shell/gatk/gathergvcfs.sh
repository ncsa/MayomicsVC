#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## gathergvcfs.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Gathers/Aggregates scattered per-sample gvcf files and creates an index 
# 
#############################################################################

 USAGE:
 gathergvcfs.sh       -s           <sample_name>
                      -b           <chr1.vcf[,chr2.vcf,...]>
                      -S           </path/to/gatk/executable>
                      -e           </path/to/java_options_file>
                      -F           </path/to/shared_functions.sh>
                      -d           turn on debug mode

 EXAMPLES:
 gathergvcfs.sh -h
 gathergvcfs.sh -s sample -b chr1.vcf,chr2.vcf,chr3.vcf -S /path/to/gatk/executable  -e /path/to/java_options_file -F /path/to/shared_functions.sh -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=gathergvcfs.sh
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
while getopts ":hs:b:S:e:F:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo -e "\n${DOCS}\n"
                        exit 0
                        ;;
                s )  # Sample name
                         SAMPLE=${OPTARG}
                         checkArg
                         ;;
                b )  # Full path to the input gvcfs or list of gvcfs
                        INPUTGVCFS=${OPTARG}
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

## Check if Sample Name variable exists
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

## Create log for JOB_ID/script
ERRLOG=${SAMPLE}.gathergvcfs.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=${SAMPLE}.gathergvcfs_gatk.log
truncate -s 0 ${TOOL_LOG}

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTGVCFS+x}" "Missing input gvcf option: -b" $LINENO
for GVCF in $(echo ${INPUTGVCFS} | sed "s/,/ /g")
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



#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------


## Defining file names
GVCFS=$( echo ${INPUTGVCFS} | sed "s/,/ --INPUT /g" | tr "\n" " " )
OUTGVCF=${SAMPLE}.g.vcf




#-------------------------------------------------------------------------------------------------------------------------------
## Gathers multiple VCF files from a scatter operation into a single VCF file. 
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Gathering gvcf variants files across a sample"

## gatk/picard MergeVcfs command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GatherVcfs aggregation error. " ' INT TERM EXIT
${GATKEXE} ${JAVA_OPTS} MergeVcfs --INPUT ${GVCFs} --OUTPUT ${OUTGVCF} >> ${TOOL_LOG} 2>&1 
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] Gathering of input GVCFs complete."




#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output gvcf and its index. Open read permissions to the user group
checkFile ${OUTGVCF} "Output variants file ${OUTGVCF} is empty." $LINENO
checkFile ${OUTGVCF}.idx "Output variants index file ${OUTGVCF}.idx is empty." $LINENO

chmod g+r ${OUTGVCF}
chmod g+r ${OUTGVCF}.idx


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
