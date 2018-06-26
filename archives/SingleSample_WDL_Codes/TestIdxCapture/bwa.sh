#!/bin/bash

REF_FA=$1
READ_1=$2
READ_2=$3
SAMPLE_NAME=$4
REF_DIR=$5

# Base name of reference.fa file
TEMP_VAR=${REF_FA}

REF_NAME=${TEMP_VAR%.*}

cd "${REF_DIR}"

if [ "${REF_NAME}*" != "${REF_FA}" ]
then
   # Create symbolic links to all reference index files
   #   
   ln -s ${REF_NAME}* ${REF_DIR}
fi

/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${SAMPLE_NAME}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${SAMPLE_NAME}" $REF_FA $READ_1 $READ_2 > /projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/IdxCaptureOutputs_temp/${SAMPLENAME}.aligned.sam

find -type l -delete
