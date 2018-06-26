#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/DesignBlock_2a/Tasks/realignment.wdl" as REALIGN
import "src/wdl_scripts/DesignBlock_2a/Tasks/bqsr.wdl" as BQSR
import "src/wdl_scripts/DesignBlock_2a/Tasks/haplotyper.wdl" as HAPLOTYPER

workflow CallBlock2aTasks {
   
   call REALIGN.realignmentTask  
   
   call BQSR.bqsrTask {
      input:
         InputAlignedSortedDedupedRealignedBam = realignmentTask.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamIdx = realignmentTask.AlignedSortedDedupedRealignedBamIdx
   
   }

   call HAPLOTYPER.variantCallingTask {
      input:
         InputAlignedSortedDedupedRealignedBam = realignmentTask.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBam = realignmentTask.AlignedSortedDedupedRealignedBamIdx,
         RecalTable = bqsrTask.RecalTable      

   }
   
}
