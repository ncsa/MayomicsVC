#!/bin/bash

Aligned_Bam=$1
sampleName=$2
Exit_Code=$3
Novosort=$4
Failure_Logs=$5

# The 'set' command is used set and examine shell options, as well as set positional parameters
set -x

# Check to see if input files are non-zero
[ -s ${Aligned_Bam} ] || echo "Input BAM File is Empty" >> ${Failure_Logs}

# Record the start of the program execution
StartTime=`date +%s`

# Novosort Tools is used to created sort BAM Files 
${NOVOSORT} -c 36 -i -o ${sampleName}.aligned.sorted.bam

# Record the End of the program execution
EndTime=`date +%s`

# Calculate the time to complete the run BWA and Samtools on each sample
echo "${sampleName} ran Novosort for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}

# The 'if' check to see if any of the samples have failed this step          
if [ $? -ne 0 ]; then
   echo '${sampleName} has failed at the Novosort Step' >> ${Exit_Code}
fi

[ ! -f ${sampleName}.aligned.sorted.bam ] && echo "aligned sorted bam not created" >> ${Failure_Logs}

