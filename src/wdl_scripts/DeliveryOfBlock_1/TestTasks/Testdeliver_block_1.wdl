#########################################################################

#       This WDL script calls the Delivery of Design Block 1 Task      ##

#########################################################################

import "src/wdl_scripts/DeliveryOfBlock_1/Tasks/deliver_block_1.wdl" as DELIVER1

workflow CallDeliverBlock1Task {

   call DELIVER1.deliverBlock1Task as DDB1

}

