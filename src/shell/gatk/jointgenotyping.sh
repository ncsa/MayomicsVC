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
# Aggregate gvcf files produced from the haplotypecaller then do joint calling
# with GenotypeGVCFs
# 
#############################################################################

 USAGE:
 jointgenotyping.sh   -b           <sample1.g.vcf[,sample2.g.vcf,...]>
                      -S           </path/to/gatk/executable>
                      -G           <reference_genome>
                      -D           <dbsnp.vcf>
                      -J           </path/to/java8_executable>
                      -e           <java_vm_options>
                      -F           </path/to/shared_functions.sh>
                      -I           <genomic_intervals>
                      -o           <extra_genotypegvcf_options>
                      -d           turn on debug mode

 EXAMPLES:
 jointgenotyping.sh -h
 jointgenotyping.sh -b sample1.g.vcf,sample2.g.vcf,sample3.g.vcf -S /path/to/gatk/executable -G reference.fa -D dbsnp.vcf -J /path/to/java8_executable -e "'-Xms2G -Xmx8G'" -F /path/to/shared_functions.sh -I chr20 -o "'--sample_ploidy 2 --useNewAFCalculator'" -d

 NOTE: In order for getops to read in a string arguments for -e (java_vm_options) or -o (extra_genotypegvcf_options), the argument needs to be quoted with a double quote (") followed by a single quote ('). See the example above.
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
while getopts ":hb:S:G:D:J:e:F:I:o:d" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
                        echo -e "\n${DOCS}\n"
                        exit 0
                        ;;
                b )  # Full path to the input gvcfs or list of gvcfs
                        INPUTGVCFS=${OPTARG}
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
                I )  # Genomic intervals overwhich to operate
                        INTERVALS=${OPTARG}
                        checkArg
                        ;;
                o ) # Extra options and arguments to GenotypeGVCF, input as a long string, can be empty if desired
                        GENOTYPEGVCF_OPTIONS=${OPTARG}
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

checkVar "${INTERVALS+x}" "Missing Intervals option: -I ${INTERVALS}" $LINENO

checkVar "${REF+x}" "Missing reference genome option: -G" $LINENO
checkFile ${REF} "Reference genome file ${REF} is empty or does not exist." $LINENO

checkVar "${DBSNP+x}" "Missing dbSNP option: -D" $LINENO
checkFile ${DBSNP} "DBSNP ${DBSNP} is empty or does not exist." $LINENO


checkVar "${GENOTYPEGVCF_OPTIONS+x}" "Missing extra GenotypeGVCFs options option: -o" $LINENO




#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME PARSING
#-------------------------------------------------------------------------------------------------------------------------------
JAVA_OPTS_PARSED=`sed -e "s/'//g" <<< ${JAVA_OPTS_STRING}`
GENOTYPEGVCF_OPTIONS_PARSED=`sed -e "s/'//g" <<< ${GENOTYPEGVCF_OPTIONS}`

## Defining file names
VARIANTS=$( echo ${INPUTGVCFS} | sed "s/,/ --variant /g" | tr "\n" " " )
OUTDB=${INTERVALS}.DB
OUTVCF=${INTERVALS}.vcf

rm -rf ${OUTDB}



#-------------------------------------------------------------------------------------------------------------------------------
## Aggregate variants across samples with GATK GenomicsDBImport
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Aggregating per sample variants across interval: ${INTERVALS}"

## gatk GenomicsDBImport command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GenomicsDBImport aggregation error. " ' INT TERM EXIT
# Memory setting based on GATK's recommendation here: https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/JointGenotypingWf.wdl#L313
${GATKEXE} --java-options "${JAVA_OPTS_PARSED}" GenomicsDBImport --genomicsdb-workspace-path ${OUTDB} --intervals ${INTERVALS} --batch-size 50 --consolidate true --variant ${VARIANTS} --reader-threads 5 -ip 500 >> ${TOOL_LOG} 2>&1 
EXITCODE=$?
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[GATKEXE] Aggregation of input GVCFs complete."




## Check for creation of output Database. 
checkDir ${OUTDB} "Output variants database ${OUTDB} is empty or malformated" $LINENO



#-------------------------------------------------------------------------------------------------------------------------------
## Perform Joint Calling with GATK GenotypeGVCFs
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[GATKEXE] Joint calling via GenotypeGVCFs"

## gatk genotypegvcf command
TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. GenotypeGVCFs joint calling error. " ' INT TERM EXIT
${GATKEXE}  --java-options "${JAVA_OPTS_PARSED}"  GenotypeGVCFs --reference ${REF} --output ${OUTVCF} --variant gendb://${OUTDB} --dbsnp ${DBSNP} ${GENOTYPEGVCF_OPTIONS_PARSED} --intervals ${INTERVALS}  >> ${TOOL_LOG} 2>&1
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
