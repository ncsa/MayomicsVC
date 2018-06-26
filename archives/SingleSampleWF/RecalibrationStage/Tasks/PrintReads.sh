#!/bin/bash

RefFasta=$1
Aligned_Sorted_Dedupped_Bam=$2
Recal_Report=$3
sampleName=$4
JAVA=$5
GATK=$6
Exit_Code=$7
Failure_Logs=$8

# Check to see if input files are non-zero
[ -s ${Aligned_Sorted_Dedupped_Bam} ] || echo "Aligned Sorted Dedupped Bam File is Empty" >> ${Failure_Logs}

${JAVA} -Xmx16g -jar ${GATK} -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Bam} -T PrintReads -BQSR ${Recal_Report} --out ${sampleName}.aligned.sorted.dedupped.realigned.recalibrated.bam -nct 6

if [ $? -ne 0 ]; then
   echo '${sampleName} has failed at the Print Reads Step' >> ${Exit_Code}
fi

[ ! -f ${sampleName}.aligned.sorted.dedupped.realigned.recalibrated.bam ] && echo "Aligned sorted dedupped realigned recalibrated bam not created" >> ${Failure_Logs}
