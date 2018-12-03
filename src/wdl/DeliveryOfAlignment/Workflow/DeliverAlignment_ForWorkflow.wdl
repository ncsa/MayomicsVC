##########################################################################################################
####              This WDL script is used to run the Delivery of alignment Block as individual module   ##
##########################################################################################################

import "src/wdl_scripts/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment

workflow CallDeliveryAlignmentTask {

   File GlobalAlignedSortedDedupedBam
   File GlobalAlignedSortedDedupedBamBai

############## BOILERPLATE FOR DELIVERY of DESIGN BLOCK 1 #######################################

   String DebugMode
   String SampleName

   File BashPreamble
   File DeliveryAlignment_Script
   String DeliveryFolder_Alignment

   File WorkflowJson

#####################################################################################          
   
   call DELIVER_Alignment.deliverAlignmentTask as deliver1 {
      input:
         SampleName = SampleName,
         AlignedSortedDedupedBam = GlobalAlignedSortedDedupedBam,
         AlignedSortedDedupedBamBai = GlobalAlignedSortedDedupedBamBai,
         WorkflowJson = WorkflowJson,
         DeliveryAlignment_Script = DeliveryAlignment_Script,
         DeliveryFolder_Alignment = DeliveryFolder_Alignment,
         DebugMode = DebugMode,
         BashPreamble = BashPreamble
   }
    
}
