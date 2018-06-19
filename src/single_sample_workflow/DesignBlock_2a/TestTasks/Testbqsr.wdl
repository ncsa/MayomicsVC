##############################################################

###       This WDL script calls the BQSR WDL Task       ##

##############################################################

import "DesignBlock_2a/Tasks/bqsr.wdl" as BQSR

workflow CallbqsrTask {

   call BQSR.bqsrTask

}
