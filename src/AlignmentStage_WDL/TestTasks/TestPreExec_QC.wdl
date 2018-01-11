###########################################################################

##          The wdl scripts check for executables and inputs          ## 

###########################################################################

import "../Tasks/PreExec_QC.wdl" as PreQC

workflow Call_PreExec_QC {

   call PreQC.PreExec_QC

}
 
