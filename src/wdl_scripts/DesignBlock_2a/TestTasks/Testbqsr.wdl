##############################################################

###       This WDL script calls the BQSR WDL Task       ##

##############################################################

import "wdl_scripts/DesignBlock_2a/Tasks/bqsr.wdl" as BQSR

workflow CallbqsrTask {

   call BQSR.bqsrTask

}
