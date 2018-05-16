#!/bin/bash

set -x

Aligned_Bam=$1
sampleName=$2
Exit_Code=$3
SORT=$4
Failure_Logs=$5 

[ -s ${Aligned_Bam} ] || echo "Input BAM File is Empty" >> ${Failure_Logs}

# Novosort Tools is used to created sort BAM Files 
${SORT} -c 36 -i -o ${sampleName}.aligned.sorted.bam ${Aligned_Bam}

# The 'if' check to see if any of the samples have failed this step          
if [ $? -ne 0 ]; then
   echo '${sampleName} has failed at the Novosort Step' >> ${Exit_Code}
fi

[ ! -f ${sampleName}.aligned.sorted.bam ] && echo "aligned sorted bam not created" >> ${Failure_Logs}
