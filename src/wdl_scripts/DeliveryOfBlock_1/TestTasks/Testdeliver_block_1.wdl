#########################################################################

#       This WDL script calls the Delivery of Design Block 1 Task      ##

#########################################################################

import "src/wdl_scripts/DeliveryOfBlock_1/Workflow/DeliverBlock1_ForWorkflow.wdl" as DWF1

workflow CallDeliverBlock1Task {

   call DWF1.DeliveryBlock1Task as DDB1

}

