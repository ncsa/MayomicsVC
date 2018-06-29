#!/bin/bash

Input_Read1=$1
Input_Read2=$2
sampleName=$3
REF_GENOME=$4
BWA=$5
SAMTOOLS=$6
Exit_Code=$7
Failure_Logs=$8

# The 'set' command is used to set and examine shell options, as well as set positional parameters
set -x

# Check to see if input files are non-zero
[[ -s ${Input_Read1} ]] || echo -e "Input Read 1 File is Empty" >> ${Failure_Logs}
[[ -s ${Input_Read2} ]] || echo -e "Input Read 2 File is Empty" >> ${Failure_Logs}

# Record the start of the program execution
StartTime=`date +%s`

# BWA Mem Tool is used to create aligned SAM file from the input FASTA File
# PIPESTATUS is an internal bash variable which holds the exit code of commands in the pipe
${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${REF_GENOME} ${Input_Read1} ${Input_Read2} | ${SAMTOOLS} view -@ 17 -bSu -> ${sampleName}.aligned.bam; B=(${dollar}{PIPESTATUS[*]})

# Record the End of the program execution
EndTime=`date +%s`

# Calculate the time to complete the run BWA and Samtools on each sample
echo -e "${sampleName} ran BWA Mem and Samtools View for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}

# Checks to see if either BWA or Samtools has a non-zero exit code
if [ ${dollar}{B[0]} -ne 0 ]
then
   echo -e "${sampleName} exited BWA with code ${dollar}{B[0]}" >> ${Exit_Code}
fi

if [ ${dollar}{B[1]} -ne 0 ]
then
   echo -e "${sampleName} exited SAMTOOLS with code ${dollar}{B[1]}" >> ${Exit_Code}
fi

# Check to see if the BAM file is created or not 
[[ ! -f ${sampleName}.aligned.bam ]] && echo -e "aligned bam not created" >> ${Failure_Logs}

