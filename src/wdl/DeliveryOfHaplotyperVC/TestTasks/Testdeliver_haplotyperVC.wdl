#########################################################################

#       This WDL script calls the Delivery of HaplotyperVC block Task     ##

#########################################################################



import "MayomicsVC/src/wdl/DeliveryOfHaplotyperVC/Tasks/deliver_HaplotyperVC.wdl" as DELIVER_HaplotyperVC

workflow CallDeliverHaplotyperVCTask {

   call DELIVER_HaplotyperVC.deliverHaplotyperVCTask as DHVCB

}

