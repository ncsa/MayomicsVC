##########################################################################################################
####           This WDL script is used to run the delivery of alignment Block as a standalone module   ###
##########################################################################################################

import "MayomicsVC/src/wdl/gatk/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment

workflow CallDeliverAlignmentTask {

   call DELIVER_Alignment.deliverAlignmentTask as DAB

}

