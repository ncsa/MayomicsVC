#!/bin/bash

set -x

Failure_Logs=$1
BWA=$2
RefFasta=$3
Input_Read1=$4
Input_Read2=$5
SAMTOOLS=$6
sampleName=$7


# Check to see if input files are non-zero
[[ -s ${Input_Read1} ]] || echo -e "Input Read 1 File is Empty" >> ${Failure_Logs}
[[ -s ${Input_Read2} ]] || echo -e "Input Read 2 File is Empty" >> ${Failure_Logs}

# Record the start of the program execution
StartTime=`date +%s`



# BWA Mem Tool is used to create aligned SAM file from the input FASTA File
# PIPESTATUS is an internal bash variable which holds the exit code of commands in the pipe

${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} | ${SAMTOOLS} view -b > ${sampleName}.aligned.bam; B=(${PIPESTATUS[*]})



# Record the End of the program execution
EndTime=`date +%s`

# Calculate the time to complete the run BWA and Samtools on each sample
echo -e "${sampleName} ran BWA MEM and samtools view for $((${EndTime}-${StartTime})) seconds" >> ${Failure_Logs}


# Checks to see if either BWA or Samtools has a non-zero exit code
if [ ${B[0]} -ne 0 ]
then
   echo -e "${sampleName} exited BWA with code ${B[0]}" >> ${Exit_Code}
fi

if [ ${B[1]} -ne 0 ]
then
   echo -e "${sampleName} exited SAMTOOLS with code ${B[1]}" >> ${Exit_Code}
fi

