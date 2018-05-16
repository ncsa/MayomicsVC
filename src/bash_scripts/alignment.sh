#!/bin/bash

################################################################################################################################
#
# Align reads using BWA-MEM. Part of the MayomicsVC Workflow.
# 
# Usage:
# trim_sequences.sh <read1.fq> <read2.fq> <reference_genome> <output_directory> </path/to/BWA> </path/to/SAMTools> <threads> </path/to/error_log>
#
################################################################################################################################

## Input and Output parameters
INPUT1=$1
INPUT2=$2
REFGEN=$3
OUTDIR=$4
BWA=$5
SAMTOOLS=$6
thr=$7
ERRLOG=$8

#set -x

## Check if input files, directories, and variables are non-zero
if [ ! -s $INPUT1 ]
then 
        echo -e "$0 stopped at line $LINENO. \nREASON=Input read 1 file $INPUT1 is empty." >> ${ERRLOG}
	exit 1;
fi
if [ ! -s $INPUT2 ]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Input read 2 file $INPUT2 is empty." >> ${ERRLOG}
	exit 1;
fi
if [ ! -s $REFGEN ]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Reference genome file $REFGEN is empty." >> ${ERRLOG}
        exit 1;
fi
if [ ! -d $OUTDIR ]
then
	echo -e "$0 stopped at line $LINENO. \nREASON=Output directory $OUTDIR does not exist." >> ${ERRLOG}
	exit 1;
fi
if [ ! -d $BWA ]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=BWA directory $BWA does not exist." >> ${ERRLOG}
	exit 1;
fi
if [ ! -d $SAMTOOLS ]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=SAMTools directory $SAMTOOLS does not exist." >> ${ERRLOG}
        exit 1;
fi
if (( $thr % 2 != 0 ))
then
	thr=$((thr-1))
fi
if [ ! -s $ERRLOG ]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Error log file $ERRLOG does not exist." >> ${ERRLOG}
        exit 1;
fi

## Parse filenames without full path
name=$(echo "$INPUT1" | sed "s/.*\///")
full=$INPUT1
sample1=${full##*/}
sample=${sample1%%.*}
OUT=${OUTDIR}/${sample}.sam
OUTBAM=${OUTDIR}/${sample}.bam
SORTBAM=${OUTDIR}/${sample}.sorted.bam

## Record start time
StartTime=`date "+%m-%d-%Y %H:%M:%S"`
echo "[BWA-MEM] START. ${StartTime}"

## BWA-MEM command, run for each read against a reference genome.
## Allocates all available threads to the process.
${BWA}/bwa mem -t $thr -M -k 32 -I 300,00 $REFGEN $INPUT1 $INPUT2 > $OUT &
wait
EndTime=`date "+%m-%d-%Y %H:%M:%S"`
echo "[BWA-MEM] Aligned reads $INPUT1 and $INPUT2 to reference $REFGEN. ${EndTime}"

## Convert SAM to BAM
echo "[SAMTools] Converting SAM to BAM..."
${SAMTOOLS}/samtools view -@ $thr -S -b $OUT > $OUTBAM &
wait
BAMEnd=`date "+%m-%d-%Y %H:%M:%S"`
echo "[SAMTools] Converted output to BAM format. ${BAMEnd}"

## Sort BAM
echo "[SAMTools] Sorting BAM..."
${SAMTOOLS}/samtools sort -@ $thr $OUTBAM -o $SORTBAM
SortEnd=`date "+%m-%d-%Y %H:%M:%S"`
echo "[SAMTools] Sorted BAM. ${SortEnd}"


## Open read permissions to the user group
chmod g+r $OUT
chmod g+r $OUTBAM
chmod g+r $SORTBAM

echo "[BWA-MEM] Finished alignment. Aligned reads found in BAM format at ${SORTBAM}. ${EndTime}"
