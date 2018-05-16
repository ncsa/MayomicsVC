#!/bin/bash

################################################################################################################################
#
# Trim input sequences using Cutadapt. Part of the MayomicsVC Workflow.
# 
# Usage:
# trim_sequences.sh <adapters.fa> <read1.fq> <read2.fq> <output_directory> </path/to/cutadapt> <threads> </path/to/error_log>
#
################################################################################################################################

## Input and Output parameters
ADAPTERS=$1
INPUT1=$2
INPUT2=$3
OUTDIR=$4
CUTADAPT=$5
thr=$6
ERRLOG=$7

#set -x

## Check if input files, directories, and variables are non-zero
if [ ! -s $ADAPTERS ]
then 
	echo -e "$0 stopped at line $LINENO. \nREASON=Adapters fasta file $ADAPTERS is empty." >> ${ERRLOG}
	exit 1;
fi
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
if [ ! -d $OUTDIR ]
then
	echo -e "$0 stopped at line $LINENO. \nREASON=Output directory $OUTDIR does not exist." >> ${ERRLOG}
	exit 1;
fi
if [ ! -d $CUTADAPT ]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Cutadapt directory $CUTADAPT does not exist." >> ${ERRLOG}
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

## Parse filename without full path
full1=$INPUT1
full2=$INPUT2
READ1=${full1##*/}
READ2=${full2##*/}
read1=${READ1%%.*}
read2=${READ2%%.*}
OUT1=${OUTDIR}/${read1}.trimmed.fq.gz
OUT2=${OUTDIR}/${read2}.trimmed.fq.gz

## Record start time
StartTime=`date "+%m-%d-%Y %H:%M:%S"`
echo "[CUTADAPT] START. ${StartTime}"

## Cutadapt command, run for each fastq and each adapter sequence in the adapter FASTA file.
## Allocates half of the available threads to each process.
${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=$((thr/2)) -o ${OUT1} $INPUT1 >> ${read1}.cutadapt.log &
${CUTADAPT}/cutadapt -a file:${ADAPTERS} --cores=$((thr/2)) -o ${OUT2} $INPUT2 >> ${read2}.cutadapt.log &
wait
echo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT logs: ${OUT}/${read1}.cutadapt.log ${OUT}/${read2}.cutadapt.log"

## Record end of the program execution
EndTime=`date "+%m-%d-%Y %H:%M:%S"`

## Open read permissions to the user group
chmod g+r $OUT1
chmod g+r $OUT2

echo "[CUTADAPT] Finished trimming adapter sequences. Trimmed reads found at ${OUT1} and ${OUT2}. ${EndTime}"
