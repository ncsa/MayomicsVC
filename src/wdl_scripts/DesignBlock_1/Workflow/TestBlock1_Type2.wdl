##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/DesignBlock_1/Tasks/trim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/DesignBlock_1/Tasks/alignment.wdl" as ALIGNMENT
import "src/wdl_scripts/DesignBlock_1/Tasks/dedup.wdl" as DEDUP 

workflow CallBlock1Tasks {
   
   call CUTADAPTTRIM.trimsequencesTask 
    
   call ALIGNMENT.alignmentTask {
      input:
         InputRead1 = trimsequencesTask.TrimmedInputRead1,
         InputRead2 = trimsequencesTask.TrimmedInputRead2
   }
   
   call DEDUP.dedupTask {
      input:
         InputAlignedSortedBam  = alignmentTask.AlignedSortedBam,
         InputAlignedSortedBamIdx = alignmentTask.AlignedSortedBamIdx
   }
    
   output {
     
      GlobalAlignedSortedDedupedBam = DEDUP.dedupTask.AlignedSortedDeduppedBam,
      GlobalAlignedSortedDedupedBamIdx = DEDUP.dedupTask.AlignedSortedDeduppedBamIdx
     
   }    

}
