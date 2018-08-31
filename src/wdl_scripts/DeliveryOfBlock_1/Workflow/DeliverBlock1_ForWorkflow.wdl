##########################################################################################################
####              This WDL script is used to run the Delivery of Design Blosk 1 as individual module    ##
##########################################################################################################

import "src/wdl_scripts/DeliveryOfBlock_1/Tasks/deliver_block_1.wdl" as DELIVER_1

workflow CallDeliveryBlock1Task {

   File GlobalAlignedSortedDedupedBam
   File GlobalAlignedSortedDedupedBamBai

############## BOILERPLATE FOR DELIVERY of DESIGN BLOCK 1 #######################################

   String DebugMode
   String SampleName

   File DeliveryBlock_1_Script
   String DeliveryFolder_Block_1


#####################################################################################          
   
   call DELIVER_1.deliverBlock1Task as deliver1 {
      input:
         SampleName = SampleName
         AlignedSortedDedupedBam = GlobalAlignedSortedDedupedBam,
         AlignedSortedDedupedBamBai = GlobalAlignedSortedDedupedBamBai,
         DeliveryBlock_1_Script = DeliveryBlock_1_Script,
         DeliveryFolder_Block_1 = DeliveryFolder_Block_1,
         DebugMode = DebugMode,
   }
    
}
