###########################################################################################
#
##            This WDL script writes out sequencing data using Print Reads               ##

##                                  Script Options                                     ##

##             -T                "Name of the tool to run"                  (Required)
###             --out            "Write output to this BAM filename"        (Required)
###             -nct             "Number of Threads"                        (Optional)
###             -R               "Reference FASTA file"                     (Optional)
############################################################################################

task PrintReadsTask {
   File RefFasta
   File Aligned_Sorted_Dedupped_Bam
   File Recal_Report
   String sampleName
   String JAVA
   String GATK
   String Exit_Code
   String Failure_Logs
   String BashScriptPath
   String dollar = "$"         

   command {
   
   /bin/bash ${BashScriptPath} ${RefFasta} ${Aligned_Sorted_Dedupped_Bam} ${Recal_Report} ${sampleName} ${JAVA} ${GATK} ${Exit_Code} ${Exit_Code}
   
   }

   output {
      File Aligned_Sorted_Dedupped_Realigned_Recalibrated_bam = "${sampeName}.aligned.sorted.dedupped.realigned.recalibrated.bam"
      String Global_sampleName = "${sampleName}"
      Array[String] Collect_Name_Bam =[Global_sampleName, Aligned_Sorted_Dedupped_Realigned_Recalibrated_bam]

   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true

   }

}
