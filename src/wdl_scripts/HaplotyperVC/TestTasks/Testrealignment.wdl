#########################################################################################################

###                This WDL script is performs Realignment on Input sorted Deduped BAM                 ##

#########################################################################################################

import "src/wdl_scripts/HaplotyperVC/Tasks/realignment.wdl" as REALIGN

workflow CallRealignmentTask {

   call REALIGN.realignmentTask

}
