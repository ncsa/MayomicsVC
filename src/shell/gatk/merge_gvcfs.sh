#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## merge_gvcfs.sh MANIFEST, USAGE DOCS, SET CHECKS
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
 merge_gvcfs.sh       -s           <sample_name>
                      -b           <chr1.vcf[,chr2.vcf,...]>
                      -S           </path/to/gatk/executable>
                      -J           </path/to/java8_executable>
                      -e           <java_vm_options>
                      -F           </path/to/shared_functions.sh>
                      -d           turn on debug mode

 EXAMPLES:
 merge_gvcfs.sh -h
 merge_gvcfs.sh -s sample -b chr1.vcf,chr2.vcf,chr3.vcf -S /path/to/gatk/executable -J /path/to/java8_executable -e "'-Xms2G -Xmx8G'" -F /path/to/shared_functions.sh -d

 NOTE: In order for getops to read in a string arguments for -e (java_vm_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=merge_gvcfs.sh
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
                b )  # Full path to the input gvcfs or list of gvcfs
                        INPUTGVCFS=${OPTARG}
                        checkArg
                        ;;
                S )  # Full path to gatk executable
                        GATKEXE=${OPTARG}
                        checkArg
                        ;;
                J ) # Path to JAVA8 exectable. The variable needs to be small letters so as not to explicitly change the user's $PATH variable
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
ERRLOG=${SAMPLE}.merge_gvcfs.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=${SAMPLE}.merge_gvcfs_gatk.log
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

## Check java8 path and options 
checkVar "${java+x}" "Missing JAVA path option: -J" $LINENO
checkFileExe ${java} "REASON=JAVA file ${java} is not executable or does not exist." $LINENO
checkVar "${JAVA_OPTS_STRING+x}" "Missing specification of JAVA memory options: -e" $LINENO



#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------

JAVA_OPTS_PARSED=`sed -e "s/'//g" <<< ${JAVA_OPTS_STRING}`

## Defining file names
GVCFS=$( echo ${INPUTGVCFS} | sed "s/,/ --INPUT /g" | tr "\n" " " )
OUTGVCF=${SAMPLE}.g.vcf




#-------------------------------------------------------------------------------------------------------------------------------
## Merge multiple gVCF files from a scatter operation into a single gVCF file. 
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Merging gvcf variants files across a sample"

## gatk/picard MergeVcfs command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. MergeVcfs aggregation error. " ' INT TERM EXIT
${GATKEXE} --java-options  "${JAVA_OPTS_PARSED}" MergeVcfs --INPUT ${GVCFS} --OUTPUT ${OUTGVCF} >> ${TOOL_LOG} 2>&1 
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
