################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/WOWScripts/WOW_Type3/BWAMemSamtoolsView.wdl" as BWA

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/WOWScripts/WOW_Type3/Novosort.wdl" as NOVOSORT

workflow WorkFlowTimes2 {

   File RefFasta

   File Ref_Amb_File
   File Ref_Dict_File
   File Ref_Ann_File
   File Ref_Bwt_File
   File Ref_Fai_File
   File Ref_Pac_File
   File Ref_Sa_File

   File InputSamplesFile

   String BWA
   String SAMTOOL

   String SORT

   call BWA.CallReadMappingTask as wf_BWA {
      input :
         RefFasta = RefFasta,
         Ref_Amb_File = Ref_Amb_File,
         Ref_Dict_File = Ref_Dict_File,
         Ref_Ann_File = Ref_Ann_File,
         Ref_Bwt_File = Ref_Bwt_File,
         Ref_Fai_File = Ref_Fai_File,
         Ref_Pac_File = Ref_Pac_File,
         Ref_Sa_File = Ref_Sa_File,
         InputSamplesFile = InputSamplesFile,
      
         BWA = BWA,
         SAMTOOL = SAMTOOL
   }


   call NOVOSORT.CallNovosortTask as wf_nsort {
      input :
         SORT = SORT,
         inputAlignedBam = wf_BWA.Global_AlignedBam
         
   }

}
