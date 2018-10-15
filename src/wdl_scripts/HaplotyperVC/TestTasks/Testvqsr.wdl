##############################################################

###       This WDL script calls the VQSR WDL Task       ##

##############################################################

import "src/wdl_scripts/HaplotyperVC/Tasks/vqsr.wdl" as BQSR

workflow CallvqsrTask {

   call VQSR.vqsrTask

}
