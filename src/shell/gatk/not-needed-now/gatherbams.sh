#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## gatherbams.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Merge bams produced in alignment 
# 
#############################################################################

 USAGE:
 gatherbams.sh     -s           <sample_name> 
                   -b           <lane1_aligned.sorted.bam[,lane2_aligned.sorted.bam,...]>
                   -S           </path/to/gatk/executable> 
                   -F           </path/to/shared_functions.sh>
                   -e           </path/to/java_options_file>
                   -d           turn on debug mode

 EXAMPLES:
 gatherbams.sh -h
 gatherbams.sh -s sample -b lane1.aligned.sorted.bam,lane2.aligned.sorted.bam,lane3.aligned.sorted.bam -S /path/to/gatk/executable -F /path/to/shared_functions.sh -e /path/to/java_options_file -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=gatherbams.sh
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
while getopts ":hs:b:S:F:e:d" OPT
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
                b )  # Full path to the input BAM or list of BAMS
                        INPUTBAM=${OPTARG}
                        checkArg
                        ;;
                S )  # Full path to gatk executive 
                        GATKEXE=${OPTARG}
                        checkArg
                        ;;
                F )  # Path to shared_functions.sh
                        SHARED_FUNCTIONS=${OPTARG}
                        checkArg
                        ;;
                 e )  # Path to file containing JAVA options to pass into the gatk command 
                        JAVA_OPTS_FILE=${OPTARG}
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
ERRLOG=${SAMPLE}.gatherbams.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.gatherbams_gatk.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTBAM+x}" "Missing input BAM option: -b" $LINENO
for LANE in $(echo ${INPUTBAM} | sed "s/,/ /g")
do
        checkFile ${LANE} "Input sorted BAM file ${LANE} is empty or does not exist." $LINENO
        checkFile ${LANE}.bai "Input sorted BAM index file ${LANE}.bai is empty or does not exist." $LINENO
done
checkVar "${GATKEXE+x}" "Missing GATKEXE path option: -S" $LINENO
checkFileExe ${GATKEXE} "REASON=SAMTOOLS file ${GATKEXE} is not an executable or does not exist." $LINENO

checkVar "${JAVA_OPTS_FILE+x}" "Missing file of JAVA path and options: -e" $LINENO
source ${JAVA_OPTS_FILE}
checkVar "${JAVA_PATH+x}" "Missing JAVA path from: -e ${JAVA_OPTS_FILE}" $LINENO
checkDir ${JAVA_PATH} "REASON=JAVA path ${JAVA_PATH} is not a directory or does not exist." $LINENO
checkVar "${JAVA_OPTS+x}" "Missing JAVA options from: -e ${JAVA_OPTS_FILE}" $LINENO




#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Defining file names
BAMS=`sed -e 's/,/ --INPUT /g' <<< ${INPUTBAM}`  ## Replace commas with spaces
MERGED_BAM=${SAMPLE}.bam
TOOL_LOG=${SAMPLE}.gatherbams_gatk.log





#-------------------------------------------------------------------------------------------------------------------------------
## BAM Merging
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Merging BAMs if list was given."

## gatk merge command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GATK BAM merging error. " ' INT TERM EXIT
${GATKEXE} ${JAVA_OPTS} GatherBamFiles --INPUT ${BAMS} --OUTPUT ${MERGED_BAM} --CREATE_INDEX true >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] BAM merging complete."





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index. Open read permissions to the user group
checkFile ${MERGED_BAM} "Output merged BAM file ${MERGED_BAM} is empty." $LINENO
checkFile ${SAMPLE}.bai "Output merged BAM index file ${SAMPLE}.bai is empty." $LINENO

chmod g+r ${MERGED_BAM}
chmod g+r ${SAMPLE}.bai



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
