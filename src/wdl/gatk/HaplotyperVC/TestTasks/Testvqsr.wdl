##############################################################

###       This WDL script calls the VQSR WDL Task       ##

##############################################################

import "MayomicsVC/src/wdl/gatk/HaplotyperVC/Tasks/vqsr.wdl" as BQSR

workflow CallvqsrTask {

   call VQSR.vqsrTask

}
