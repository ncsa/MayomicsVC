##############################################################

###       This WDL script calls the VQSR WDL Task       ##

##############################################################

import "MayomicsVC/src/wdl/HaplotyperVC/Tasks/vqsr.wdl" as BQSR

workflow CallvqsrTask {

   call VQSR.vqsrTask

}
