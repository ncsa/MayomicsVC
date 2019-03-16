#########################################################################

##       This WDL script calls the Bam Merging  WDL Task     ##

#########################################################################

import "MayomicsVC/src/wdl/gatk/Alignment/Tasks/merge_aligned_bam.wdl" as MERGE 

workflow Callmerge_aligned_bamTask {

   call MERGE.mergebamTask

}

