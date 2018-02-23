################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/WOWScripts/BWAMemSamtoolsView.wdl" as BWA

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/WOWScripts/Novosort.wdl" as NOVOSORT

workflow wf_Master {
   
   call BWA.CallReadMappingTask as wf_BWA 

   call NOVOSORT.CallNovosortTask {
      input :
         inputAlignedBam = wf_BWA.Global_AlignedBam
         
   }
}
