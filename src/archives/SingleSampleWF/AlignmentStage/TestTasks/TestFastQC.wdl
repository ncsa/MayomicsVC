#######################################################################################################

##            This WDL script performs Quality Control on input FastQ files            ##

#######################################################################################################

import "AlignmentStage_WDL/Tasks/FastQC.wdl" as FASTQC

workflow CallFastqQualityControlTask {

      call FASTQC.FastqQualityControlTask


} # End of Workflow block
