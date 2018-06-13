#########################################################################################################

###                This WDL script is performs Realignment on Input sorted Deduped BAM                 ##

#########################################################################################################

import "DesignBlock_1/Tasks/Sentieon/realignment.wdl" as REALIGN

workflow CallRealignmentTask {

   call REALIGN.realignmentTask

}
