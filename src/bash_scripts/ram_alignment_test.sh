#!/bin/bash

################################################################################################################################
#
# Align reads using BWA-MEM. Part of the MayomicsVC Workflow.
# 
# Usage:
# alignment.sh -s <sample_name> -r <read1.fq> -R <read2.fq> -G <reference_genome> -O <output_directory> -B </path/to/BWA> -T </path/to/SAMTools> -t <threads> -SE <Is_single_end?> -e </path/to/error_log>
#
################################################################################################################################

## Input and Output parameters
while getopts ":h:s:r:R:G:O:B:T:t:SE:e:" OPT
do
        case ${OPT} in
                h )
                        echo "Usage:"
                        echo "  bash alignment.sh -h       Display this help message."
                        echo "  bash alignment.sh [-s <sample_name>] [-r <read1.fq>] [-R <read2.fq>] [-G <reference_genome>] [-O <output_directory>] [-B </path/to/BWA>] [-T </path/to/SAMTools>] [-t threads] [-SE single-end? (true/false)] [-e </path/to/error_log>] "
                        ;;
                s )
                        s=${OPTARG}
                        echo $s
                        ;;
                r )
                        r=${OPTARG}
                        echo $r
                        ;;
                R )
                        R=${OPTARG}
                        echo $R
                        ;;
                G )
                        G=${OPTARG}
                        echo $G
                        ;;
                O )
                        O=${OPTARG}
                        echo $O
                        ;;
                B )
                        B=${OPTARG}
                        echo $B
                        ;;
		T )
			T=${OPTARG}
			echo $T
			;;
                t )
                        t=${OPTARG}
                        echo $t
                        ;;
                SE )
                        SE=${OPTARG}
                        echo $SE
                        ;;
                e )
                        e=${OPTARG}
                        echo $e
                        ;;
        esac
done



INPUT1=${r}
INPUT2=${R}
SAMPLE=${s}
REFGEN=${G}
OUTDIR=${O}
BWA=${B}
SAMTOOLS=${T}
THR=${t}
IS_SINGLE_END=${SE}
ERRLOG=${e}

#set -x

## Check if input files, directories, and variables are non-zero
if [[ ! -s ${INPUT1} ]]
then 
        echo -e "$0 stopped at line $LINENO. \nREASON=Input read 1 file ${INPUT1} is empty." >> ${ERRLOG}
	exit 1;
fi
if [[ ${IS_SINGLE_END} == false ]]
then
        if [[ ! -s ${INPUT2} ]]
        then
                echo -e "$0 stopped at line $LINENO. \nREASON=Input read 2 file ${INPUT2} is empty." >> ${ERRLOG}
	        exit 1;
        fi
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
if [[ ! -f ${ERRLOG} ]]
then
        echo -e "$0 stopped at line $LINENO. \nREASON=Error log file ${ERRLOG} does not exist." >> ${ERRLOG}
        exit 1;
fi

## Parse filenames without full path
name=$(echo "${INPUT1}" | sed "s/.*\///")
full=${INPUT1}
sample1=${full##*/}
sample=${sample1%%.*}
OUT=${OUTDIR}/${SAMPLE}.sam
OUTBAM=${OUTDIR}/${SAMPLE}.bam
SORTBAM=${OUTDIR}/${SAMPLE}.sorted.bam

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
