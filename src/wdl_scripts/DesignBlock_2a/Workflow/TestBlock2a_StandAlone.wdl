#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/DesignBlock_2a/Tasks/realignment.wdl" as REALIGNMENT
import "src/wdl_scripts/DesignBlock_2a/Tasks/bqsr.wdl" as BQSR
import "src/wdl_scripts/DesignBlock_2a/Tasks/haplotyper.wdl" as HAPLOTYPER
import "src/wdl_scripts/DesignBlock_2a/Tasks/vqsr.wdl" as VQSR

workflow CallBlock2aTasks {
   
   call REALIGNMENT.realignmentTask  as realign 
   
   call BQSR.bqsrTask as bqsr {
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamIdx = realign.AlignedSortedDedupedRealignedBamIdx,
   }

   call HAPLOTYPER.variantCallingTask as haplotype { 
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBamIdx,
   }

   call VQSR.vqsrTask as vqsr {
      input:
         InputVCF = haplotype.VCF,
         InputVCFIdx = haplotype.VcfIdx,
   }

   
}
