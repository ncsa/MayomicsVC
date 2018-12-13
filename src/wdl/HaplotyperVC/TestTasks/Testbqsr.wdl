##############################################################

###       This WDL script calls the BQSR WDL Task       ##

##############################################################

import "MayomicsVC/src/wdl/HaplotyperVC/Tasks/bqsr.wdl" as BQSR

workflow CallbqsrTask {

   call BQSR.bqsrTask

}
