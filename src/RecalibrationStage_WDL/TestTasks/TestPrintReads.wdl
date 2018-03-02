#########################################################################################################

###  This WDL script performs Print Reads to write out sequence read data   ###

#########################################################################################################

import "RecalibrationStage_WDL/Tasks/PrintReads.wdl" as PRINTREADS

workflow CallPrintReadsTask {

   File RefFasta
   File Aligned_Sorted_Dedupped_Bam
   File Recal_Report
   String sampleName
   String JAVA
   String GATK
   String Exit_Code
   String Failure_Logs

   # The 2-D array stores information of all the samples 
   Array[Array[File]] inputRecalReportDeduppedBam

   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(Bam in inputRecalReportDeduppedBam) {

      # BWA Mem is included as a sub task and it is called inside the workflow
      call PRINTREADS.PrintReadsTask {
         input :
            sampleName = sample[0],
            Recal_Report = sample[1],
            Aligned_Sorted_Dedupped_Bam = sample[2]
      }

   } # End of scatter block

   output {
      Array[File] Global_Aligned_Sorted_Dedupped_Realigned_Recalibrated_bam = PrintReadsTask.Collect_Name_Bam
   }

} # End of Workflow block

