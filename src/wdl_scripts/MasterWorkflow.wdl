###############################################################################################
####              This WDL script is used to run Alignment and HaplotyperVC blocks together  ##
###############################################################################################

import "src/wdl_scripts/Alignment/Workflow/Alignment_ForWorkflow.wdl" as WF1
import "src/wdl_scripts/DeliveryOfAlignment/Workflow/DeliverAlignment_ForWorkflow.wdl" as DWF1
import "src/wdl_scripts/HaplotyperVC/Workflow/HaplotyperVC_ForWorkflow.wdl" as WF2a
import "src/wdl_scripts/DeliveryOfHaplotyperVC/Workflow/DeliverHaplotyperVC_ForWorkflow.wdl" as DWF2a

workflow MasterWF {

   call WF1.CallAlignmentTasks as AlignmentBlock

   call DWF1.CallDeliveryAlignmentTask as DAB {
      input:
         GlobalAlignedSortedDedupedBam = AlignmentBlock.GlobalAlignedSortedDedupedBam,
         GlobalAlignedSortedDedupedBamBai = AlignmentBlock.GlobalAlignedSortedDedupedBamBai
   }

   call WF2a.CallHaplotyperVCTasks as HaplotyperVCBlock {
      input:
         GlobalAlignedSortedDedupedBam = AlignmentBlock.GlobalAlignedSortedDedupedBam,
         GlobalAlignedSortedDedupedBamBai = AlignmentBlock.GlobalAlignedSortedDedupedBamBai
   }

   call DWF2a.CallDeliveryHaplotyperVCTask as DHVCB {
      input:
         GlobalRecalibratedVcf = HaplotyperVCBlock.GlobalRecalibratedVcf,
         GlobalRecalibratedVcfIdx = HaplotyperVCBlock.GlobalRecalibratedVcfIdx
   }
}
