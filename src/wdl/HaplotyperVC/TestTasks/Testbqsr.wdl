##############################################################

###       This WDL script calls the BQSR WDL Task       ##

##############################################################

import "src/wdl/HaplotyperVC/Tasks/bqsr.wdl" as BQSR

workflow CallbqsrTask {

   call BQSR.bqsrTask

}
