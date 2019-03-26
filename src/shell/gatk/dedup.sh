#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## dedup.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Deduplicate BAMs with Picard. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 dedup.sh          -s           <sample_name> 
                   -b           <aligned_sorted_merged.bam>
                   -S           </path/to/gatk/executable> 
                   -J           </path/to/java8_executable>
                   -e           <java_vm_options>
                   -F           </path/to/shared_functions.sh>
                   -d           turn on debug mode

 EXAMPLES:
 dedup.sh -h
 dedup.sh -s sample -b aligned_sorted_merged.bam -S /path/to/gatk/executable -J /path/to/java8_executable -e "'-Xms2G -Xmx8G'" -F /path/to/shared_functions.sh -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=dedup.sh
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
while getopts ":hs:b:S:J:e:F:d" OPT
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
                b )  # Full path to the input BAM
                        INPUTBAM=${OPTARG}
                        checkArg
                        ;;
                S )  # Full path to gatk executable
                        GATKEXE=${OPTARG}
                        checkArg
                        ;;
                J ) # Path to JAVA8 exectable. The variable needs to be small letters so as not to explicitly change the user's $PATH variabl
                        java=${OPTARG}
                        checkArg
                        ;;
                e ) # JAVA options string to pass into the gatk command 
                        JAVA_OPTS_STRING=${OPTARG}
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
ERRLOG=${SAMPLE}.dedup.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.dedup_picard.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"


## Check java8 path and options 
checkVar "${java+x}" "Missing JAVA path option: -J" $LINENO
checkFileExe ${java} "REASON=JAVA file ${java} is not executable or does not exist." $LINENO
checkVar "${JAVA_OPTS_STRING+x}" "Missing specification of JAVA memory options: -e" $LINENO


## Check if input files, directories, and variables are non-zero
checkVar "${INPUTBAM+x}" "Missing input BAM option: -b" $LINENO
checkFile ${INPUTBAM} "Input sorted BAM file ${INPUTBAM} is empty or does not exist." $LINENO
checkFile ${INPUTBAM}.bai "Input sorted BAM index file ${INPUTBAM}.bai is empty or does not exist." $LINENO

checkVar "${GATKEXE+x}" "Missing GATK path option: -S" $LINENO
checkFileExe ${GATKEXE} "REASON=GATK file ${GATKEXE} is not executable or does not exist." $LINENO




#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Defining file names
OUT=${SAMPLE}.bam
DEDUPMETRICS=${SAMPLE}.dedup_metrics.txt
TOOL_LOG=${SAMPLE}.dedup_picard.log

JAVA_OPTS_PARSED=`sed -e "s/'//g" <<< ${JAVA_OPTS_STRING}`




#-------------------------------------------------------------------------------------------------------------------------------
## DEDUPLICATION
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[PICARD] Deduplicating BAM."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Picard Deduplication error. " ' INT TERM EXIT
${GATKEXE} --java-options "${JAVA_OPTS_PARSED}" MarkDuplicates --INPUT ${INPUTBAM} --METRICS_FILE ${DEDUPMETRICS} --OUTPUT ${OUT} >> ${TOOL_LOG}  2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[PICARD] Deduplication Finished. Deduplicated BAM found at ${OUT}"




#-------------------------------------------------------------------------------------------------------------------------------
## BAM INDEXONG 
#-------------------------------------------------------------------------------------------------------------------------------

## Index BAM 
logInfo "[PICARD] Indexing BAM..."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Picard BAM indexing error. " ' INT TERM EXIT
${GATKEXE} --java-options  "${JAVA_OPTS_PARSED}" BuildBamIndex --INPUT ${OUT} --OUTPUT ${OUT}.bai >> ${TOOL_LOG} 2>&1
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[PICARD] Indexed BAM output."




#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index. Open read permissions to the user group
checkFile ${OUT} "Output deduplicated BAM file ${OUT} is empty." $LINENO
checkFile ${OUT}.bai "Output deduplicated BAM index file ${OUT}.bai is empty." $LINENO

chmod g+r ${OUT}
chmod g+r ${OUT}.bai
chmod g+r ${DEDUPMETRICS}




#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
