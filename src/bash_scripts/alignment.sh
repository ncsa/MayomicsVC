#!/bin/bash

################################################################################################################################
#
# Align reads using BWA-MEM. Part of the MayomicsVC Workflow.
# 
# Usage:
# trim_sequences.sh <read1.fq> <read2.fq> <reference_genome> <output_directory> </path/to/BWA> </path/to/SAMTools> <threads> <Is_single_end?> </path/to/error_log>
#
################################################################################################################################

## Input and Output parameters
INPUT1=$1
INPUT2=$2
SAMPLE=$3
REFGEN=$4
OUTDIR=$5
BWA=$6
SAMTOOLS=$7
THR=$8
IS_SINGLE_END=$9
ERRLOG=$10

#set -x

## Check if input files, directories, and variables are non-zero
if [[ ! -s ${INPUT1} ]]
then 
        echo -e "$0 stopped at line $LINENO. \nREASON=Input read 1 file ${INPUT1} is empty." >> ${ERRLOG}
	exit 1;
fi
if [[ ! -s ${INPUT2} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Input read 2 file ${INPUT2} is empty." >> ${ERRLOG}
	exit 1;
fi
if [[ ! -s ${REFGEN} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Reference genome file ${REFGEN} is empty." >> ${ERRLOG}
        exit 1;
fi
if [[ ! -d ${OUTDIR} ]]
then
	echo -e "$0 stopped at line $LINENO. \nREASON=Output directory ${OUTDIR} does not exist." >> ${ERRLOG}
	exit 1;
fi
if [[ ! -d ${BWA} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=BWA directory ${BWA} does not exist." >> ${ERRLOG}
	exit 1;
fi
if [[ ! -d ${SAMTOOLS} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=SAMTools directory ${SAMTOOLS} does not exist." >> ${ERRLOG}
        exit 1;
fi
if (( ${THR} % 2 != 0 ))
then
	THR=$((THR-1))
fi
if [[ ! -s ${ERRLOG} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Error log file ${ERRLOG} does not exist." >> ${ERRLOG}
        exit 1;
fi

## Parse filenames without full path
name=$(echo "${INPUT1}" | sed "s/.*\///")
full=${INPUT1}
sample1=${full##*/}
sample=${sample1%%.*}
OUT=${OUTDIR}/${sample}.sam
OUTBAM=${OUTDIR}/${sample}.bam
SORTBAM=${OUTDIR}/${sample}.sorted.bam

## Record start time
START_TIME=`date "+%m-%d-%Y %H:%M:%S"`
echo "[BWA-MEM] START. ${START_TIME}"

## BWA-MEM command, run for each read against a reference genome.
## Allocates all available threads to the process.
if [[ ${IS_SINGLE_END} == true ]]
then
	${BWA}/bwa mem -t ${THR} -M -k 32 ${REFGEN} ${INPUT1} > ${OUT} &
	wait
else
	${BWA}/bwa mem -t ${THR} -M -k 32 -I 300,30 ${REFGEN} ${INPUT1} ${INPUT2} > ${OUT} &
	wait
fi
END_TIME=`date "+%m-%d-%Y %H:%M:%S"`
echo "[BWA-MEM] Aligned reads ${SAMPLE} to reference ${REFGEN}. ${END_TIME}"

## Convert SAM to BAM
echo "[SAMTools] Converting SAM to BAM..."
${SAMTOOLS}/samtools view -@ ${THR} -S -b ${OUT} > ${OUTBAM} &
wait
BAM_TIME=`date "+%m-%d-%Y %H:%M:%S"`
echo "[SAMTools] Converted output to BAM format. ${BAM_TIME}"

## Sort BAM
echo "[SAMTools] Sorting BAM..."
${SAMTOOLS}/samtools sort -@ ${THR} ${OUTBAM} -o ${SORTBAM}
SORT_TIME=`date "+%m-%d-%Y %H:%M:%S"`
echo "[SAMTools] Sorted BAM. ${SORT_TIME}"


## Open read permissions to the user group
chmod g+r ${OUT}
chmod g+r ${OUTBAM}
chmod g+r ${SORTBAM}

echo "[BWA-MEM] Finished alignment. Aligned reads found in BAM format at ${SORTBAM}. ${END_TIME}"
