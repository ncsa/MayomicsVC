#######################################################################################

##       This WDL script is performs Novosort to sort the input BAM Files     ##

#######################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/MultiSample_WDL_Codes/Tasks/Novosort.wdl" as NSORT

workflow Call_Novosort {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   #The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)
   
   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {

      # Novosort is included as a sub task and it is called inside the workflow
      call NSORT.Novosort {
         input :
            sampleName = sample[0]
      }
  }
}
   
