##########################################################################################################
####              This WDL script is used to run the Delivery of the HaplotyperVC module    ##
##########################################################################################################

import "src/wdl_scripts/DeliveryOfHaplotyperVC/Tasks/deliver_HaplotyperVC.wdl" as DELIVER_HaplotyperVC

workflow CallDeliveryHaplotyperVCTask {

   File GlobalRecalibratedVcf
   File GlobalRecalibratedVcfIdx


############## BOILERPLATE FOR DELIVERY of HaplotyperVC #######################################

   String DebugMode
   String SampleName

   File DeliveryHaplotyperVC_Script
   String DeliveryFolder_HaplotyperVC

   File WorkflowJson
   File BashPreamble

#####################################################################################          
   
   call DELIVER_HaplotyperVC.deliverHaplotyperVCTask as deliverHaplotyperVC {
      input:
         BashPreamble = BashPreamble,
         SampleName = SampleName,
         RecalibratedVcf = GlobalRecalibratedVcf,
         RecalibratedVcfIdx = GlobalRecalibratedVcfIdx,
         WorkflowJson = WorkflowJson,
         DeliveryHaplotyperVC_Script = DeliveryHaplotyperVC_Script,
         DeliveryFolder_HaplotyperVC = DeliveryFolder_HaplotyperVC,
         DebugMode = DebugMode
   }
    
}
