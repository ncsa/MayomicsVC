##########################################################################################################
####              This WDL script is used to run the Delivery of Design Blosk 2a as individual module    ##
##########################################################################################################

import "src/wdl_scripts/DeliveryOfBlock_2a/Tasks/deliver_block_2a.wdl" as DELIVER_2a

workflow CallDeliveryBlock2aTask {

   File GlobalRecalibratedSnpVcf
   File GlobalRecalibratedSnpVcfIdx
   File GlobalRecalibratedIndelVcf
   File GlobalRecalibratedIndelVcfIdx


############## BOILERPLATE FOR DELIVERY of DESIGN BLOCK 2a #######################################

   String DebugMode

   File DeliveryBlock_2a_Script
   String DeliveryFolder_Block_2a

   File WorkflowJson

#####################################################################################          
   
   call DELIVER_2a.deliverBlock2aTask as deliver2a {
      input:
         RecalibratedSnpVcf = GlobalRecalibratedSnpVcf,
         RecalibratedSnpVcfIdx = GlobalRecalibratedSnpVcfIdx,
         RecalibratedIndelVcf = GlobalRecalibratedIndelVcf,
         RecalibratedIndelVcfIdx = GlobalRecalibratedIndelVcfIdx,
         WorkflowJson = WorkflowJson,
         DeliveryBlock_2a_Script = DeliveryBlock_2a_Script,
         DeliveryFolder_Block_2a = DeliveryFolder_Block_2a,
         DebugMode = DebugMode
   }
    
}
