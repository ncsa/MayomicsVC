#########################################################################################################
###           This WDL script is used to run the Recalibration steps as individual modules             ##
#########################################################################################################

import "RecalibrationStage_WDL/TestTasks/TestBaseRecalibrator.wdl" as BASERECALIBRATOR
import "RecalibrationStage_WDL/TestTasks/TestPrintReads.wdl" as PRINTREADS


workflow CallRecalibrationStageTasks {

   call BASERECALIBRATOR.CallBaseRecalibrationTask as wf_BaseRecalibrate

   call PRINTREADS.CallPrintReadsTask {
      input :
         inputRecalReportDeduppedBam = wf_BaseRecalibrate.Global_Recal_Report
   }

}

