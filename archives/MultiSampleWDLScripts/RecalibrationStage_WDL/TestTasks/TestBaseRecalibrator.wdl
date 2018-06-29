#########################################################################################################

###  This WDL script performs Base Recalibration to detects systematic errors in base quality scores   ##

#########################################################################################################

import "RecalibrationStage_WDL/Tasks/BaseRecalibrator.wdl" as BASERECALIBRATION

workflow CallBaseRecalibrationTask {
   
   File RefFasta
   File Aligned_Sorted_Dedupped_Bam
   File Millsand1000GIndels
   String JAVA
   String GATK
   String Exit_Code
   String Failure_Logs

   # The 2-D array stores information of all the Aligned Sorted Dedupped Bam
   Array[Array[File]] inputAlignedSortedDeduppedBam

   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(Bam in inputAlignedSortedDeduppedBam) {

      # BWA Mem is included as a sub task and it is called inside the workflow
      call BASERECALIBRATION.BaseRecalibrationTask {
         input :
            sampleName = Bam[0],
            Aligned_Sorted_Dedupped_Bam = Bam[1]
            
      }

   } # End of scatter block

   output {
      Array[File] Global_Recal_Report = BaseRecalibrationTask.Collect_Name_Report
   }
} # End of Workflow block

