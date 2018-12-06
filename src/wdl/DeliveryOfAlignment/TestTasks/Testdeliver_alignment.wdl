#########################################################################

#       This WDL script calls the Delivery of alignment block Task     ##

#########################################################################

import "src/wdl/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment

workflow CallDeliverAlignmentTask {

   call DELIVER_Alignment.deliverAlignmentTask as DAB

}

