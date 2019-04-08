#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "MayomicsVC/src/wdl/sentieon/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "MayomicsVC/src/wdl/sentieon/HaplotyperVC/Tasks/bqsr.wdl" as BQSR
import "MayomicsVC/src/wdl/sentieon/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER
import "MayomicsVC/src/wdl/sentieon/HaplotyperVC/Tasks/vqsr.wdl" as VQSR

workflow CallHaplotyperVCTasks {
   
   call REALIGNMENT.realignmentTask  as realign 
   
   call BQSR.bqsrTask as bqsr {
      input:
         InputBams = realign.OutputBams,
         InputBais = realign.OutputBais,
   }

   call HAPLOTYPER.variantCallingTask as haplotype { 
      input:
         InputBams = realign.OutputBams,
         InputBais = realign.OutputBais,
         RecalTable = bqsr.RecalTable,
   }

   call VQSR.vqsrTask as vqsr {
      input:
         InputVcf = haplotype.OutputVcf,
         InputVcfIdx = haplotype.OutputVcfIdx,
   }

   
}
