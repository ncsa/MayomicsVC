#########################################################################################################
###           This WDL script is used to run the Recalibration steps as individual modules             ##
#########################################################################################################

import "RecalibrationStage_WDL/TestTasks/TestBaseRecalibrator.wdl" as BASERECALIBRATOR
import "RecalibrationStage_WDL/TestTasks/TestPrintReads.wdl" as PRINTREADS


workflow CallRecalibrationStageTasks {

   call BASERECALIBRATOR.CallBaseRecalibrationTask

   call PRINTREADS.CallPrintReadsTask 

}

