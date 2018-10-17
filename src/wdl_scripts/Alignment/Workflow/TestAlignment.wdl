##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/Alignment/Tasks/trim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/Alignment/Tasks/alignment.wdl" as ALIGNMENT
import "src/wdl_scripts/Alignment/Tasks/dedup.wdl" as DEDUP 

workflow CallAlignmentTasks {
   
   call CUTADAPTTRIM.trimsequencesTask as trimseq 
    
   call ALIGNMENT.alignmentTask as align {
      input:
         InputRead1 = trimseq.TrimmedInputRead1,
         InputRead2 = trimseq.TrimmedInputRead2
   }
   
   call DEDUP.dedupTask as dedup {
      input:
         InputAlignedSortedBam  = align.AlignedSortedBam,
         InputAlignedSortedBamBai = align.AlignedSortedBamBai
   }
    
}
