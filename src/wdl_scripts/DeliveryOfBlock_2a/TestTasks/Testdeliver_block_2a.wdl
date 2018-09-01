#########################################################################

#       This WDL script calls the Delivery of Design Block 2a Task     ##

#########################################################################



import "src/wdl_scripts/DeliveryOfBlock_2a/Tasks/deliver_block_2a.wdl" as DELIVER2a

workflow CallDeliverBlock2aTask {

   call DELIVER2a.deliverBlock1Task as DDB2a

}

