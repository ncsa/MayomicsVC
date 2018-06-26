#########################################################################################################

###                This WDL script is performs Realignment on Input sorted Deduped BAM                 ##

#########################################################################################################

import "src/wdl_scripts/DesignBlock_2a/Tasks/realignment.wdl" as REALIGN

workflow CallRealignmentTask {

   call REALIGN.realignmentTask

}
