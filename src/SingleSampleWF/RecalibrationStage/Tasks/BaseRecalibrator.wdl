###########################################################################################

##                     This WDL script performs Base Recalibration                       ##

##                                  Script Options                                     ##

##             -T               "Name of the tool to run"                  (Required)
##             --out            "output recalibration table file"          (Required)
##             -nct             "Number of Threads"                        (Optional)
##             -R               "Reference FASTA file"                     (Optional)
##             -knownsites      "A database of known polymorphic sites     (Optional)
###########################################################################################

# The Task block is where the variables and the functions are defined for performing a certain task

task BaseRecalibrationTask {

   File RefFasta
   File Aligned_Sorted_Dedupped_Bam
   File Millsand1000GIndels
   String sampleName
   String JAVA
   String GATK
   String Exit_Code
   String Failure_Logs
   String BashScriptPath
   String dollar = "$"

   command {

      /bin/bash ${BashScriptPath} ${RefFasta} ${Aligned_Sorted_Dedupped_Bam} ${Millsand1000GIndels} ${sampleName} ${JAVA} ${GATK} ${Exit_Code} ${Failure_Logs}

   }

   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs   

   output {
      File Recal_Report = "${sampleName}_recal_report.grp"
      String Global_sampleName = "${sampleName}"
      File Global_Aligned_Sorted_Dedupped_Bam = "${Aligned_Sorted_Dedupped_Bam}"
      Array[String] Collect_Name_Report = [Global_sampleName, Recal_Report, Global_Aligned_Sorted_Dedupped_Bam]
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true

   }

}
