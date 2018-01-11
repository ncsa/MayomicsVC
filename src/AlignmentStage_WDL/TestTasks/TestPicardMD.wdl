###############################################################################################

#       This WDL script is performs PICARD Mark Duplicates on input sorted BAM Files     ##

###############################################################################################

import "../Tasks/PicardMD.wdl" as PICARD

workflow Call_PicardMD {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   #The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {

      # Novosort is included as a sub task and it is called inside the workflow
      call PICARD.Picard_MarkDuplicates {
         input :
            sampleName = sample[0]
      }
   }
}
