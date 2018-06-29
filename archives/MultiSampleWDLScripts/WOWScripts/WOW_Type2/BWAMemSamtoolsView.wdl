###########################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/MayomicsVC/archives/MultiSampleWDLScripts/WOWScripts/WOW_Type2/BWA.wdl" as BWA 

workflow CallReadMappingTask {

   File RefFasta

   File Ref_Amb_File
   File Ref_Dict_File
   File Ref_Ann_File
   File Ref_Bwt_File
   File Ref_Fai_File
   File Ref_Pac_File
   File Ref_Sa_File

   String BWA                 
   String SAMTOOL      

   File InputSamplesFile

   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   scatter(sample in inputsamples) {
   
      call BWA.ReadMappingTask as B {
         input :
            sampleName = sample[0],
            Input_Read1 = sample[1],
            Input_Read2 = sample[2],

            RefFasta = RefFasta,
            Ref_Amb_File = Ref_Amb_File,
            Ref_Dict_File = Ref_Dict_File,
            Ref_Ann_File = Ref_Ann_File,
            Ref_Bwt_File = Ref_Bwt_File,
            Ref_Fai_File = Ref_Fai_File,
            Ref_Pac_File = Ref_Pac_File,
            Ref_Sa_File = Ref_Sa_File,

            BWA = BWA,
            SAMTOOL = SAMTOOL 
      }
   
   }

   output {
      Array[Array[String]] Global_AlignedBam = B.Collect_AlignedBam
   }
   
}
