##############################################################

###       This WDL script calls the VQSR WDL Task       ##

##############################################################

import "src/wdl_scripts/DesignBlock_2a/Tasks/vqsr.wdl" as BQSR

workflow CallvqsrTask {

   call VQSR.vqsrTask

}
