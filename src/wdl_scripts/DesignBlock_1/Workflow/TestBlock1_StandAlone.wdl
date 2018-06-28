##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/DesignBlock_1/Tasks/trim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/DesignBlock_1/Tasks/alignment.wdl" as ALIGNMENT
import "src/wdl_scripts/DesignBlock_1/Tasks/dedup.wdl" as DEDUP 

workflow CallBlock1Tasks {
   
   call CUTADAPTTRIM.trimsequencesTask as trimseq {
      input:
         InputRead1 = InputRead1,
         InputRead2 = InputRead2
   }
    
   call ALIGNMENT.alignmentTask as align {
      input:
         InputRead1 = trimseq.TrimmedInputRead1,
         InputRead2 = trimseq.TrimmedInputRead2
   }
   
   call DEDUP.dedupTask as dedup {
      input:
         InputAlignedSortedBam  = align.AlignedSortedBam,
         InputAlignedSortedBamIdx = align.AlignedSortedBamIdx
   }
    
}
