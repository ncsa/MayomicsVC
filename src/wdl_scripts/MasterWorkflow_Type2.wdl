###############################################################################################
####              This WDL script is used to run Design Block 1 and 2a together              ##
###############################################################################################

import "src/wdl_scripts/DesignBlock_1/Workflow/TestBlock1_Type2.wdl" as WF1
import "src/wdl_scripts/DesignBlock_2a/Workflow/TestBlock2a_Type2.wdl" as WF2

workflow MasterWF {

   call WF1.CallBlock1Tasks as CALLWF1

   call WF2.CallBlock2aTasks as CALLWF2 {
      input:
         GlobalAlignedSortedDedupedBam = CALLWF1.GlobalAlignedSortedDedupedBam,
         GlobalAlignedSortedDedupedBamIdx = CALLWF1.GlobalAlignedSortedDedupedBamIdx

   }

}
