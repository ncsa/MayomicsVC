###########################################################################

##          The wdl scripts check for executables and inputs          ## 

###########################################################################

import "AlignmentStage_WDL/Tasks/PreExecQC.wdl" as PREEXECQC

workflow CallQualityControlTask {

   call PREEXECQC.QualityControlTask

}
 
