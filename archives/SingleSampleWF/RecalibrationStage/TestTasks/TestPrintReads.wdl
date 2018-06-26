#########################################################################################################

###  This WDL script performs Print Reads to write out sequence read data   ###

#########################################################################################################

import "RecalibrationStage_WDL/Tasks/PrintReads.wdl" as PRINTREADS

workflow CallPrintReadsTask {

      call PRINTREADS.PrintReadsTask 

} # End of Workflow block

