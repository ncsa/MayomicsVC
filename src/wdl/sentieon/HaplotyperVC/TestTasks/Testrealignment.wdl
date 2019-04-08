#########################################################################################################

###                This WDL script is performs Realignment on Input sorted Deduped BAM                 ##

#########################################################################################################

import "MayomicsVC/src/wdl/sentieon/HaplotyperVC/Tasks/realignment.wdl" as REALIGN

workflow CallRealignmentTask {

   call REALIGN.realignmentTask

}
