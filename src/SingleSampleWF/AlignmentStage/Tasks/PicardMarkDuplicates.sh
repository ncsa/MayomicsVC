#!/bin/bash

Aligned_Sorted_Bam=$1
sampleName=$2
Exit_Code=$3
JAVA=$4
PICARD=$5
Failure_Logs=$6

# The 'set' command is used to set and exaine shell options, as well as positional parameters
set -x

# Check to see if input files are non-zero
[ -s ${Aligned_Sorted_Bam} ] || echo "Aligned Sorted Bam File is Empty" >> ${Failure_Logs}

# Record the start of the program execution
StartTime=`date +%s`

# Picard Mark Duplicates is used to mark duplicates on input sorted BAMs
${JAVA} -Xmx2g -jar ${PICARD} MarkDuplicates I=${Aligned_Sorted_Bam} O=${sampleName}.aligned.sorted.dedupped.bam M=${sampleName}.PicardMetrics ASSUME_SORTED=true CREATE_INDEX=true

EndTime=`date +%s`

# Calculate the time to complete the run BWA and Samtools on each sample
echo "${sampleName} ran Picard Mark Duplicates for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}

if [ $? -ne 0 ]; then
   echo '${sampleName} has failed at the Mark Duplicates Step' >> ${Exit_Code}
fi

[ ! -f ${sampleName}.aligned.sorted.dedupped.bam ] && echo "aligned sorted dedupped bam not created" >> ${Failure_Logs}
