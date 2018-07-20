###############################################################################################
####              This WDL script is used to run Design Block 1 and 2a together              ##
###############################################################################################

import "src/wdl_scripts/DesignBlock_1/Workflow/TestBlock1.wdl" as WF1
import "src/wdl_scripts/DesignBlock_2a/Workflow/TestBlock2a.wdl" as WF2

workflow MasterWF {

   call WF1.CallBlock1Tasks as DB1

   call WF2.CallBlock2aTasks as DB2a {
      input:
         GlobalAlignedSortedDedupedBam = DB1.GlobalAlignedSortedDedupedBam,
         GlobalAlignedSortedDedupedBamIdx = DB1.GlobalAlignedSortedDedupedBamIdx
   }

}
