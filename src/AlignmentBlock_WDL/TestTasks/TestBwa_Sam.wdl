#########################################################################################################

###       This WDL script is performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/Bwa_Sam.wdl" as BWA

workflow Call_BWA {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile


   # The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

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

   } # End of scatter block

} # End of Workflow block

