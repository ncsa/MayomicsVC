#!/bin/bash

################################################################################################################################
#
# Deduplicate BAM using Sentieon Locus Collector and Dedup algorithms. Part of the MayomicsVC Workflow.
# 
# Usage:
# dedup.sh -s <sample_name> -b <aligned.sorted.bam> -O <output_directory> -S </path/to/sentieon> -t <threads> -e </path/to/error_log>
#
################################################################################################################################

## Input and Output parameters
while getopts ":h:s:b:O:S:t:e:d:" OPT
do
        case ${OPT} in
                h )  # Flag to display usage 
			echo " "
                        echo "Usage:"
			echo " "
                        echo "  bash dedup.sh -h       Display this help message."
                        echo "  bash dedup.sh [-s sample_name] [-b <aligned.sorted.bam>] [-O <output_directory>] [-S </path/to/sentieon>] [-t threads] [-e </path/to/error_log>] "
			echo " "
			exit 0;
                        ;;
		s )  # Sample name. String variable invoked with -s
			SAMPLE=${OPTARG}
			echo ${SAMPLE}
			;;
                b )  # Full path to the input BAM file. String variable invoked with -b
                        INPUTBAM=${OPTARG}
                        echo ${INPUTBAM}
                        ;;
                O )  # Output directory. String variable invoked with -O
                        OUTDIR=${OPTARG}
                        echo ${OUTDIR}
                        ;;
                S )  # Full path to sentieon directory. String variable invoked with -S
                        SENTIEON=${OPTARG}
                        echo ${SENTIEON}
                        ;;
                t )  # Number of threads available. Integer invoked with -t
                        THR=${OPTARG}
                        echo ${THR}
                        ;;
                e )  # Full path to error log file. String variable invoked with -e
                        ERRLOG=${OPTARG}
                        echo ${ERRLOG}
                        ;;
                d )  # Turn on debug mode. Boolean variable [true/false] which initiates 'set -x' to print all text
                        DEBUG=${OPTARG}
                        echo ${DEBUG}
                        ;;
        esac
done

## Turn on Debug Mode to print all code
if [[ ${DEBUG} == true ]]
then
        set -x
fi

SCRIPT_NAME=dedup.sh

## Check if input files, directories, and variables are non-zero
if [[ ! -s ${INPUTBAM} ]]
then 
        echo -e "$0 stopped at line $LINENO. \nREASON=Input sorted BAM file ${INPUTBAM} is empty." >> ${ERRLOG}
	exit 1;
fi
if [[ ! -d ${OUTDIR} ]]
then
	echo -e "$0 stopped at line $LINENO. \nREASON=Output directory ${OUTDIR} does not exist." >> ${ERRLOG}
	exit 1;
fi
if [[ ! -d ${SENTIEON} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Sentieon directory ${SENTIEON} does not exist." >> ${ERRLOG}
	exit 1;
fi
if (( ${THR} % 2 != 0 ))
then
	THR=$((THR-1))
fi
if [[ ! -f ${ERRLOG} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Error log file ${ERRLOG} does not exist." >> ${ERRLOG}
        exit 1;
fi

samplename=${SAMPLE}
SCORETXT=${OUTDIR}/${SAMPLE}.score.txt
OUT=${OUTDIR}/${SAMPLE}.deduped.bam
OUTBAMIDX=${OUTDIR}/${SAMPLE}.deduped.bam.bai
DEDUPMETRICS=${OUTDIR}/${SAMPLE}.dedup_metrics.txt

## Record start time
START_TIME=`date "+%m-%d-%Y %H:%M:%S"`
echo "[SENTIEON] Collecting info to deduplicate BAM with Locus Collector. ${START_TIME}"

## Locus Collector command
${SENTIEON}/bin/sentieon driver -t ${THR} -i ${INPUTBAM} --algo LocusCollector --fun score_info ${SCORETXT} 
wait
echo "[SENTIEON] Locus Collector finished; starting Dedup."

## Dedup command (Note: optional --rmdup flag will remove duplicates; without, duplicates are marked but not removed)
${SENTIEON}/bin/sentieon driver -t ${THR} -i ${INPUTBAM} --algo Dedup --score_info ${SCORETXT} --metrics ${DEDUPMETRICS} ${OUT}

## Record end of the program execution
END_TIME=`date "+%m-%d-%Y %H:%M:%S"`
echo "[SENTIEON] Deduplication Finished. ${END_TIME}"
echo "[SENTIEON] Deduplicated BAM found at ${OUT}"

## Open read permissions to the user group
chmod g+r ${OUT}
chmod g+r ${OUTBAMIDX}
chmod g+r ${DEDUPMETRICS}
