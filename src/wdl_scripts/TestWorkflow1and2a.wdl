###############################################################################################
####              This WDL script is used to run Design Block 1 and 2a together              ##
###############################################################################################

import "src/wdl_scripts/DesignBlock_1/Workflow/DesignBlock1_ForWorkflow.wdl" as WF1
import "src/wdl_scripts/DeliveryOfBlock_1/Workflow/DeliverBlock1_ForWorkflow.wdl" as DWF1
import "src/wdl_scripts/DesignBlock_2a/Workflow/DesignBlock2a_ForWorkflow.wdl" as WF2a
import "src/wdl_scripts/DeliveryOfBlock_2a/Workflow/DeliverBlock2a_ForWorkflow.wdl" as DWF2a

workflow MasterWF {

   call WF1.CallBlock1Tasks as DB1

   call DWF1.CallDeliveryBlock1Task as DDB1 {
      input:
         GlobalAlignedSortedDedupedBam = DB1.GlobalAlignedSortedDedupedBam,
         GlobalAlignedSortedDedupedBamBai = DB1.GlobalAlignedSortedDedupedBamBai
   }

   call WF2a.CallBlock2aTasks as DB2a {
      input:
         GlobalAlignedSortedDedupedBam = DB1.GlobalAlignedSortedDedupedBam,
         GlobalAlignedSortedDedupedBamBai = DB1.GlobalAlignedSortedDedupedBamBai
   }

   call DWF2a.CallDeliveryBlock2aTask as DDB2a {
      input:
         GlobalRecalibratedVcf = DB2a.GlobalRecalibratedVcf,
         GlobalRecalibratedVcfIdx = DB2a.GlobalRecalibratedVcfIdx,
   }
}
