#!/bin/bash

Input_Read1=$1
Input_Read2=$2
sampleName=$3
FastQCDir=$4
FASTQC=$5
JAVA=$6
Failure_Logs=$7
Exit_Code=$8

# Check to see if input files are non-zero
[ -s ${Input_Read1} ] || echo "Input Read1 FastQ is Empty" >> ${Failure_Logs}
[ -s ${Input_Read2} ] || echo "Input Read2 FastQ is Empty" >> ${Failure_Logs}

#FastQC takes a FastQ file and runs a series of tests on it to generate a comprehensive QC report
${FASTQC} --extract -j ${JAVA} -o ${FastQCDir} ${Input_Read1} ${Input_Read2}

if [ $? -ne 0 ]; then
   echo "${sampleName} has failed at the FASTQ/BAM File Quality Control Step" >> ${Exit_Code}
fi

[ ! -d ${FastQCDir} ] && echo "FASTQC directory has not been created" >> ${Failure_Logs}
