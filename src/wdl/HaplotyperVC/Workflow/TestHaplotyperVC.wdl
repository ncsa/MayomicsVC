#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "src/wdl_scripts/HaplotyperVC/Tasks/bqsr.wdl" as BQSR
import "src/wdl_scripts/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER
import "src/wdl_scripts/HaplotyperVC/Tasks/vqsr.wdl" as VQSR

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
