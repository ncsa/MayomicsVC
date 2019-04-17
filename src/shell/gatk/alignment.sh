#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------------------------------
## alignment.sh MANIFEST, USAGE DOCS, SET CHECKS
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
# Align sequences using BWA-MEM. Part of the MayomicsVC Workflow.
# 
#############################################################################

 USAGE:
 alignment.sh      -s           <sample_name> 
                   -p	    	<platform>
                   -L           <library>
                   -f           <flowcell_ID/platform_unit>
                   -c           <sequencing_center>
                   -P           paired-end reads (true/false)
                   -l           <read1.fq> 
                   -r           <read2.fq>
                   -G           <reference_genome>
                   -e           </path/to/bwa/executable>
                   -K           <chunk_size_in_bases> 
                   -o           <additional_bwa_options>
                   -S           </path/to/samtools/executable>
                   -t           <threads> 
                   -F           </path/to/shared_functions.sh>
                   -d           turn on debug mode

 EXAMPLES:
 alignment.sh -h
 alignment.sh -s sample -p platform -L library -f flowcell_ID -c center_name -l read1.fq -r read2.fq -G reference.fa -K 10000000 -o "'-M'" -S /path/to/samtools/executable -t 12 -P true -e /path/to/bwa/executable -F /path/to/shared_functions.sh -d

 NOTES: To prevent different results due to thread count, set -K to 10000000 as recommended by the functional equivalence guidelines (Regier et al 2018).
        In order for getops to read in a string arguments for -o (additional_bwa_options), the argument needs to be quoted with a double quote (") followed by a single quote (').

#############################################################################

DOCS







set -o errexit
set -o pipefail
set -o nounset

SCRIPT_NAME=alignment.sh
SGE_JOB_ID=TBD   # placeholder until we parse job ID
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
while getopts ":hs:g:p:L:f:c:l:r:G:K:o:S:t:P:e:F:d" OPT
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
                p )  # Sequencing platform
                        PLATFORM=${OPTARG}
            			checkArg
                        ;;
        		L )  # Library name
		            	LIBRARY=${OPTARG}
            			checkArg
			            ;;
        		f )  # Platform unit / flowcell ID
		            	PLATFORM_UNIT=${OPTARG}
            			checkArg
			            ;;
        		c )  # Sequencing center name
		            	CENTER_NAME=${OPTARG}
            			checkArg
			            ;;
                l )  # Full path to input read 1
                        INPUT1=${OPTARG}
            			checkArg
                        ;;
                r )  # Full path to input read 2
                        INPUT2=${OPTARG}
			            checkArg
                        ;;
                G )  # Full path to referance genome fasta file
                        REFGEN=${OPTARG}
            			checkArg
                        ;;
        		K )  # Chunk size in bases (10000000 to prevent different results based on thread count)
		            	CHUNK_SIZE=${OPTARG}
            			checkArg
			            ;;
        		o )  # Additional BWA MEM options
		            	BWA_OPTS=${OPTARG}
            			checkArg
			            ;;
        		S )  # Full path to SAMTOOLSEXE
                        SAMTOOLSEXE=${OPTARG}
		            	checkArg
                        ;;
                t )  # Number of threads available
                        THR=${OPTARG}
            			checkArg
                        ;;
                P )  # Is this a paired-end process? [true/false] Invoked with -P
                        IS_PAIRED_END=${OPTARG}
			            checkArg
                        ;;
        		e )  # Full path to BWAEXE
                        BWAEXE=${OPTARG}
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
ERRLOG=${SAMPLE}.alignment.${SGE_JOB_ID}.log
truncate -s 0 "${ERRLOG}"
truncate -s 0 ${SAMPLE}.align_bwa.log

## Write manifest to log
echo "${MANIFEST}" >> "${ERRLOG}"

## Check if input files, directories, and variables are non-zero
checkVar "${INPUT1+x}" "Missing read 1 option: -l" $LINENO
checkFile ${INPUT1} "Input read 1 file ${INPUT1} is empty or does not exist." $LINENO
checkVar "${INPUT2+x}" "Missing read 2 option: -r. If running a single-end job, set -r null in command." $LINENO

checkVar "${IS_PAIRED_END+x}" "Missing paired-end option: -P" $LINENO

if [[ "${IS_PAIRED_END}" != true ]] && [[ "${IS_PAIRED_END}" != false ]]
then
	EXITCODE=1
	logError "$0 stopped at line ${LINENO}. \nREASON=Incorrect argument for paired-end option -P. Must be set to true or false."
fi
if [[ "${IS_PAIRED_END}" == true ]]
then
        # checkFile ${INPUT2} "Input read 2 file ${INPUT2} is empty or does not exist. If running a single-end job, set -r null in command." $LINENO
	if [[ "${INPUT2}" == null ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}/ \nREASON=User specified Paired End option -P, but set read 2 option -r to null."
	fi
fi
if [[ "${IS_PAIRED_END}" == false ]]
then
	if [[  "${INPUT2}" != null ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}/ \nREASON=User specified Single End option, but did not set read 2 option -r to null."
	fi
fi

checkVar "${REFGEN+x}" "Missing reference genome option: -G" $LINENO
checkFile ${REFGEN} "Reference genome file ${REFGEN} is empty or does not exist." $LINENO

checkVar "${CHUNK_SIZE+x}" "Missing read group option: -K\nSet -K 10000000 to prevent different results based on thread count." $LINENO
checkVarInt "${CHUNK_SIZE}" "Not integer value for chunk size option: -K\nSet -K 10000000 to prevent different results based on thread count." $LINENO
if [[ ${CHUNK_SIZE} != 10000000 ]]
then
	logWarn "[BWA-MEM] Chunk size option -K set to ${CHUNK_SIZE}. When this option is not set to 10000000, there may be different results per run based on different thread counts."
fi

checkVar "${LIBRARY+x}" "Missing sequencing library option: -L" $LINENO 
checkVar "${PLATFORM+x}" "Missing platform/sequencing technology option: -p" $LINENO
checkVar "${PLATFORM_UNIT+x}" "Missing platform unit / flowcell ID option: -f" $LINENO
checkVar "${CENTER_NAME+x}" "Missing sequencing center name option: -c" $LINENO
checkVar "${BWA_OPTS+x}" "Missing additional BWA MEM options option: -O" $LINENO
checkVar "${BWAEXE+x}" "Missing BWAEXE path option: -e" $LINENO
checkFileExe ${BWAEXE} "REASON=BWA file ${BWAEXE} is not an executable or does not exist." $LINENO
checkVar "${SAMTOOLSEXE+x}" "Missing SAMTOOLSEXE path option: -S" $LINENO
checkFileExe ${SAMTOOLSEXE} "REASON=SAMTOOLS file ${SAMTOOLSEXE} is not an executable or does not exist." $LINENO
checkVar "${THR+x}" "Missing threads option: -t" $LINENO
checkVarInt "${THR}" "Not integer value for number of threads: -t" $LINENO




#-------------------------------------------------------------------------------------------------------------------------------
## FILENAME/OPTIONS PARSING
#-------------------------------------------------------------------------------------------------------------------------------

## Set output file names
OUTSAM=${SAMPLE}.sam
OUTBAM=${SAMPLE}-unsorted.bam
SORTBAM=${SAMPLE}.bam
SORTBAMIDX=${SAMPLE}.bam.bai
TOOL_LOG=${SAMPLE}.align_bwa.log

## Parse extra options if specified
BWA_OPTS_PARSED=`sed -e "s/'//g" <<< ${BWA_OPTS}`

## Parse read group name
if [[ "${IS_PAIRED_END}" == false ]]
then
	GROUP=$(basename ${INPUT1} .fastq.gz)  ## Removes path and hard-coded file extension
else
	GROUP=$(basename ${INPUT1} .fastq.gz)_$(basename ${INPUT2} .fastq.gz)
fi

#-------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------
## READ ALIGNMENT
#-------------------------------------------------------------------------------------------------------------------------------

## Record start time
logInfo "[BWA-MEM] START."

## BWA-MEM command, run for each read against a reference genome. 
if [[ "${IS_PAIRED_END}" == false ]] # Align single read to reference genome
then
	TRAP_LINE=$(($LINENO + 1))
	trap 'logError " $0 stopped at line ${TRAP_LINE}. BWA-MEM error in read alignment. " ' INT TERM EXIT
	${BWAEXE} mem ${BWA_OPTS_PARSED} -Y -R "@RG\tID:${GROUP}\tPU:${PLATFORM_UNIT}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIBRARY}\tCN:${CENTER_NAME}" -K ${CHUNK_SIZE} -t ${THR} ${REFGEN} ${INPUT1} > ${OUTSAM} # 2>>${TOOL_LOG}
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

        checkExitcode ${EXITCODE} $LINENO
else 
        # Paired-end reads aligned
	TRAP_LINE=$(($LINENO + 1))
	trap 'logError " $0 stopped at line ${TRAP_LINE}. BWA-MEM error in read alignment. " ' INT TERM EXIT
	${BWAEXE} mem ${BWA_OPTS_PARSED} -Y -R "@RG\tID:$GROUP\tPU:${PLATFORM_UNIT}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIBRARY}\tCN:${CENTER_NAME}" -K ${CHUNK_SIZE} -t ${THR} ${REFGEN} ${INPUT1} ${INPUT2} > ${OUTSAM} # 2>>${TOOL_LOG} 
	EXITCODE=$?  # Capture exit code
	trap - INT TERM EXIT

        checkExitcode ${EXITCODE} $LINENO
fi

checkFile ${OUTSAM} "Output SAM ${OUTSAM} is empty." $LINENO
logInfo "[BWA-MEM] Aligned reads ${SAMPLE} to reference ${REFGEN}."






#-------------------------------------------------------------------------------------------------------------------------------
## BAM CONVERSION
#-------------------------------------------------------------------------------------------------------------------------------

## Convert SAM to BAM
logInfo "[SAMTOOLS] Converting SAM to BAM..."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Samtools BAM conversion error. " ' INT TERM EXIT
${SAMTOOLSEXE} view -@ ${THR} -bS ${OUTSAM} -o ${OUTBAM} >> ${TOOL_LOG} 2>&1
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[SAMTOOLS] Converted output to BAM format."





#-------------------------------------------------------------------------------------------------------------------------------
## BAM SORTING
#-------------------------------------------------------------------------------------------------------------------------------

## Sort BAM
logInfo "[SAMTOOLS] Sorting BAM..."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Samtools BAM sorting error. " ' INT TERM EXIT
${SAMTOOLSEXE} sort -@ ${THR} -o ${SORTBAM} ${OUTBAM} >> ${TOOL_LOG} 2>&1
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[SAMTOOLS] Sorted output BAM."


#-------------------------------------------------------------------------------------------------------------------------------
## BAM INDEXONG 
#-------------------------------------------------------------------------------------------------------------------------------

## Index BAM 
logInfo "[SAMTOOLS] Indexing BAM..."

TRAP_LINE=$(($LINENO + 1))
trap 'logError " $0 stopped at line ${TRAP_LINE}. Samtools BAM indexing error. " ' INT TERM EXIT
${SAMTOOLSEXE} index -b ${SORTBAM} >> ${TOOL_LOG} 2>&1
EXITCODE=$?  # Capture exit code
trap - INT TERM EXIT

checkExitcode ${EXITCODE} $LINENO
logInfo "[SAMTOOLS] Indexed BAM output."





#-------------------------------------------------------------------------------------------------------------------------------
## POST-PROCESSING
#-------------------------------------------------------------------------------------------------------------------------------

## Check if BAM and index were created. Open read permissions to the user group
checkFile ${SORTBAM} "Output sorted BAM ${SORTBAM} is empty." $LINENO
checkFile ${SORTBAMIDX} "Output index for the sorted BAM ${SORTBAMIDX} is empty." $LINENO

chmod g+r ${OUTSAM}
chmod g+r ${OUTBAM}
chmod g+r ${SORTBAM}
chmod g+r ${SORTBAMIDX}

logInfo "[BWA-MEM] Finished alignment. Aligned reads found in BAM format at ${SORTBAM}."

rm ${OUTSAM} ${OUTBAM}

#-------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
## END
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
exit 0;
