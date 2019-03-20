#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## merge_bams.sh MANIFEST, USAGE DOCS, SET CHECKS
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
 merge_bams.sh     -s           <sample_name> 
                   -b           <lane1_aligned.sorted.bam[,lane2_aligned.sorted.bam,...]>
                   -S           </path/to/samtools/executable> 
                   -F           </path/to/shared_functions.sh>
                   -d           turn on debug mode

 EXAMPLES:
 merge_bams.sh -h
 merge_bams.sh -s sample -b lane1.aligned.sorted.bam,lane2.aligned.sorted.bam,lane3.aligned.sorted.bam -S /path/to/samtools/executable -F /path/to/shared_functions.sh -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=merge_bams.sh
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
while getopts ":hs:b:S:F:d" OPT
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
                S )  # Full path to samtools executive 
                        SAMTOOLSEXE=${OPTARG}
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
ERRLOG=${SAMPLE}.merge_bams.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.merge_bams_samtools.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTBAM+x}" "Missing input BAM option: -b" $LINENO
for LANE in $(echo ${INPUTBAM} | sed "s/,/ /g")
do
        checkFile ${LANE} "Input sorted BAM file ${LANE} is empty or does not exist." $LINENO
        checkFile ${LANE}.bai "Input sorted BAM index file ${LANE}.bai is empty or does not exist." $LINENO
done
checkVar "${SAMTOOLSEXE+x}" "Missing SAMTOOLSEXE path option: -S" $LINENO
checkFileExe ${SAMTOOLSEXE} "REASON=SAMTOOLS file ${SAMTOOLSEXE} is not an executable or does not exist." $LINENO





#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Defining file names
BAMS=`sed -e 's/,/ /g' <<< ${INPUTBAM}`  ## Replace commas with spaces
MERGED_BAM=${SAMPLE}.bam
TOOL_LOG=${SAMPLE}.merge_bams_samtools.log





#-------------------------------------------------------------------------------------------------------------------------------
## BAM Merging
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[SAMTOOLSEXE] Merging BAMs if list was given."

## samtools merge command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. SAMTOOLS BAM merging error. " ' INT TERM EXIT
${SAMTOOLSEXE} merge -f ${MERGED_BAM} ${BAMS} >> ${TOOL_LOG} 2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[SAMTOOLSEXE] BAM merging complete."





#-------------------------------------------------------------------------------------------------------------------------------
## BAM INDEXONG 
#-------------------------------------------------------------------------------------------------------------------------------

## Index BAM 
logInfo "[SAMTOOLS] Indexing BAM..."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Picard BAM indexing error. " ' INT TERM EXIT
${SAMTOOLSEXE} index -b ${MERGED_BAM} >> ${TOOL_LOG} 2>&1
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[SAMTOOLS] Indexed BAM output."






#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output BAM and index. Open read permissions to the user group
checkFile ${MERGED_BAM} "Output merged BAM file ${MERGED_BAM} is empty." $LINENO
checkFile ${MERGED_BAM}.bai "Output merged BAM index file ${MERGED_BAM} is empty." $LINENO

chmod g+r ${MERGED_BAM}
chmod g+r ${MERGED_BAM}.bai



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
