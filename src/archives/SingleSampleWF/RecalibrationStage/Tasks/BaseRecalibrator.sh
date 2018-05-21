#!/bin/bash

RefFasta=$1
Aligned_Sorted_Dedupped_Bam=$2
Millsand1000GIndels=$3
sampleName=$4
JAVA=$5
GATK=$6
Exit_Code=$7
Failure_Logs=$8

# Check to see if input files are non-zero
[ -s ${Aligned_Sorted_Dedupped_Bam} ] || echo "Aligned Sorted Dedupped Bam File is Empty" >> ${Failure_Logs}

# Base Recalibration detects systematic errors in base quality scores  
${JAVA} -Xmx16g -jar $GATK -T BaseRecalibrator -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Bam} -knownSites ${Millsand1000GIndels} --out ${sampleName}_recal_report.grp -nct 17

if [ $? -ne 0 ]; then
   echo '${sampleName} has failed at the Base Recalibration Step' >> ${Exit_Code}
fi

[ ! -f ${sampleName}_recal_report.grp ] && echo "Real Report not created" >> ${Failure_Logs}
