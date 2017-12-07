#########################################################################################################

###              This WDL script is used to run the Alignment steps as individual modules              ##

#########################################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/MultiSample_WDL_Codes/Tasks/Bwa_Sam.wdl" as BWA

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/MultiSample_WDL_Codes/Tasks/Novosort.wdl" as NSORT

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/MultiSample_WDL_Codes/Tasks/PicardMD.wdl" as PICARD

workflow AlignBlock_Run {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   # The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   String Capture_Exit_Code
   
   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {
   
      # BWA Mem is included as a sub task and it is called inside the workflow
      call BWA.BWA_Mem {
         input :
            sampleName = sample[0],         
            Input_Read1 = sample[1],         
            Input_Read2 = sample[2]
      }   

      call NSORT.Novosort {
         input :
            sampleName = sample[0],
            Aligned_Bam = BWA_Mem.Aligned_Bam,
            Exit_Code = Capture_Exit_Code
      }

      call PICARD.Picard_MarkDuplicates {
         input :
            sampleName = sample[0],
            Aligned_Sorted_Bam = Novosort.Aligned_Sorted_Bam,
            Exit_Code = Capture_Exit_Code
      }      

   } # End of scatter block
  
} # End of Workflow block



