#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## realignment.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Realign reads using GATK Realigner. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 realignment.sh    -s           <sample_name> 
                   -b           <sorted.deduped.bam>
                   -G           <reference_genome>
                   -k           <known_sites> (omni.vcf, hapmap.vcf, indels.vcf, dbSNP.vcf) 
                   -S           </path/to/gatk/executable> 
                   -t           <threads> 
                   -e           </path/to/file/containing/java_options>
                   -F           </path/to/shared_functions.sh>
                   -d           turn on debug mode

 EXAMPLES:
 realignment.sh -h
 realignment.sh -s sample -b sorted.deduped.bam -G reference.fa -k known1.vcf,known2.vcf,...knownN.vcf -S /path/to/gatk/executable -t 12 -e /path/to/file/containing/java/options -F /path/to/shared_functions.sh -d

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=realignment.sh
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
while getopts ":hs:b:G:k:S:t:e:F:d" OPT
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
        		b )  # Full path to the input deduped BAM
                        INPUTBAM=${OPTARG}
                        checkArg
                        ;;
                G )  # Full path to reference genome fasta file
                        REFGEN=${OPTARG}
                        checkArg
                        ;;
        		k )  # Full path to known sites files
                        KNOWN=${OPTARG}
                        checkArg
                        ;;
                S )  # Full path to gatk executable 
                        GATKEXE=${OPTARG}
                        checkArg
                        ;;
                t )  # Number of threads available
                        THR=${OPTARG}
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
ERRLOG=${SAMPLE}.realignment.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.realign_gatk.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## source the file with java path and options variables
checkVar "${JAVA_OPTS_FILE+x}" "Missing file of JAVA path and options: -e" $LINENO
source ${JAVA_OPTS_FILE}
checkVar "${JAVA_PATH+x}" "Missing JAVA path from: -e ${JAVA_OPTS_FILE}" $LINENO
checkDir ${JAVA_PATH} "REASON=JAVA path ${JAVA_PATH} is not a directory or does not exist." $LINENO
checkVar "${JAVA_OPTS+x}" "Missing JAVA options from: -e ${JAVA_OPTS_FILE}" $LINENO

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTBAM+x}" "Missing input BAM option: -b" $LINENO
checkFile ${INPUTBAM} "Input sorted BAM file ${INPUTBAM} is empty or does not exist." $LINENO
checkFile ${INPUTBAM}.bai "Input sorted BAM index file ${INPUTBAM}.bai is empty or does not exist." $LINENO

checkVar "${REFGEN+x}" "Missing reference genome option: -G" $LINENO
checkFile ${REFGEN} "Reference genome file ${REFGEN} is empty or does not exist." $LINENO

checkVar "${KNOWN+x}" "Missing known sites option ${KNOWN}: -k" $LINENO
checkVar "${GATKEXE+x}" "Missing GATK path option: -S" $LINENO
checkFileExe ${GATKEXE} "REASON=GATK file ${GATKEXE} is not executable or does not exist." $LINENO
checkVar "${THR+x}" "Missing threads option: -t" $LINENO





#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME AND OPTION PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Parse known sites list of multiple files. Create multiple -k flags for gatk
SPLITKNOWN=`sed -e 's/,/ -k /g' <<< ${KNOWN}`
echo ${SPLITKNOWN}

## Parse filenames without full path
OUT=${SAMPLE}.bam
TOOL_LOG=${SAMPLE}.





#-------------------------------------------------------------------------------------------------------------------------------
## REALIGNMENT STAGE
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[Realigner] START. Realigning deduped BAM. Using known sites at ${KNOWN} ."

## GATK Realigner command.
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GATK Realignment error. " ' INT TERM EXIT
${GATKEXE}/bin/gatk driver -t ${THR} -r ${REFGEN} -i ${INPUTBAM} --algo Realigner -k ${SPLITKNOWN} ${OUT} >> ${SAMPLE}.realign_gatk.log 2>&1
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[Realigner] Realigned reads ${SAMPLE} to reference ${REFGEN}. Realigned BAM located at ${OUT}."






#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of realigned BAM and index. Open read permissions to the user group
checkFile ${OUT} "Realigned BAM file ${OUT} is empty." $LINENO
checkFile ${OUT}.bai "Realigned BAM index file ${OUT}.bai is empty." $LINENO

chmod g+r ${OUT}
chmod g+r ${OUT}.bai






#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
