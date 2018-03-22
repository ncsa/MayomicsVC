#!/bin/bash

Input_Read1=$1
Input_Read2=$2
sampleName=$3
RefFasta=$4
Ref_Amb_File=$5
Ref_Dict_File=$6
Ref_Ann_File=$7
Ref_Bwt_File=$9
Ref_Fai_File=$10
Ref_Pac_File=$11
Ref_Sa_File=$12
BWA=$13

   set -x

   [[ -s ${Input_Read1} ]] || echo "Input Read 1 File is Empty" >> ${Failure_Logs}
   [[ -s ${Input_Read2} ]] || echo "Input Read 2 File is Empty" >> ${Failure_Logs}

   echo "Hi" 

   ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2}
   
   echo "bye"
