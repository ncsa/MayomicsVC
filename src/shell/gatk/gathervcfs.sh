#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## gathervcfs.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Gathers/Aggregates vcf files (from Picard tools)
# 
#############################################################################

 USAGE:
 gathervcfs.sh        -b           <chr1.vcf[,chr2.vcf,...]>
                      -S           </path/to/gatk/executable>
                      -F           </path/to/shared_functions.sh>
                      -J           </path/to/java8_executable>
                      -e           <java_vm_options>
                      -d           turn on debug mode

 EXAMPLES:
 gathervcfs.sh -h
 gathervcfs.sh -b chr1.vcf,chr2.vcf,chr3.vcf -S /path/to/gatk/executable -F /path/to/shared_functions.sh -J /path/to/java8_executable -e "'-Xms2G -Xmx8G'" -d

 NOTE: In order for getops to read in a string arguments for -e (java_vm_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.

#############################################################################

DOCS






set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=gathervcfs.sh
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
while getopts ":hb:S:J:e:F:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo -e "\n${DOCS}\n"
                        exit 0
                        ;;
                b )  # Full path to the input gvcfs or list of gvcfs
                        INPUTVCFS=${OPTARG}
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

## Create log for JOB_ID/script
ERRLOG=vcfs.gather.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
TOOL_LOG=vcfs.gather_gatk.log
truncate -s 0 ${TOOL_LOG}

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUTVCFS+x}" "Missing input gvcf option: -b" $LINENO
for VCF in $(echo ${INPUTVCFS} | sed "s/,/ /g")
do
        checkFile ${VCF} "Input variants file ${VCF} is empty or does not exist." $LINENO
        checkFile ${VCF}.idx "Input variants index file ${VCF}.idx is empty or does not exist." $LINENO
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
VCFS=$( echo ${INPUTVCFS} | sed "s/,/ --INPUT /g" | tr "\n" " " )
OUTVCF=GenomicGermlineVariants.vcf




#-------------------------------------------------------------------------------------------------------------------------------
## Gathers multiple VCF files from a scatter operation into a single VCF file. 
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Gathering the per interval variants files across samples"

## gatk/picard GatherVcfs command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GatherVcfs aggregation error. " ' INT TERM EXIT
${GATKEXE} --java-options "${JAVA_OPTS_PARSED}" GatherVcfs --INPUT ${VCFS} --OUTPUT ${OUTVCF} >> ${TOOL_LOG} 2>&1 
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] Gathering of input VCFs complete."




#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check for creation of output vcf. Open read permissions to the user group
checkFile ${OUTVCF} "Output variants file ${OUTVCF} is empty." $LINENO

chmod g+r ${OUTVCF}


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
