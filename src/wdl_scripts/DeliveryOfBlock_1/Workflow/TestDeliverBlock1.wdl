##########################################################################################################
####           This WDL script is used to run the delivery of Design Block 1 as a standalone module    ###
##########################################################################################################

import "src/wdl_scripts/DeliveryOfBlock_1/Tasks/deliver_block_1.wdl" as DELIVER1

workflow CallDeliverBlock1Task {

   call DELIVER1.deliverBlock1Task as DDB1

}

