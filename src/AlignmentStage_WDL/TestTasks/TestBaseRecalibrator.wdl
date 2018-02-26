#########################################################################################################

###  This WDL script performs Base Recalibration to detects systematic errors in base quality scores   ##

#########################################################################################################

import "AlignmentStage_WDL/Tasks/BaseRecalibrator.wdl" as BASERECALIBRATION

workflow CallBaseRecalibration {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   # The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {

      # BWA Mem is included as a sub task and it is called inside the workflow
      call BASERECALIBRATION.BaseRecalibration {
         input :
            sampleName = sample[0],
            Aligned_Sorted_Dedupped_Realigned_Bam = sample[1]
            
      }

   } # End of scatter block

} # End of Workflow block

