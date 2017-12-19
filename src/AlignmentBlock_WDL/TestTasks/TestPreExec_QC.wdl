###########################################################################

##          The wdl scripts check for executables and inputs          ## 

###########################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/PreExec_QC.wdl" as PreQC

workflow Call_PreExec_QC {

   call PreQC.PreExec_QC

}
 
