#########################################################################################################

###  This WDL script performs Base Recalibration to detects systematic errors in base quality scores   ##

#########################################################################################################

import "RecalibrationStage_WDL/Tasks/BaseRecalibrator.wdl" as BASERECALIBRATION

workflow CallBaseRecalibrationTask {
   
   call BASERECALIBRATION.BaseRecalibrationTask

} # End of Workflow block

