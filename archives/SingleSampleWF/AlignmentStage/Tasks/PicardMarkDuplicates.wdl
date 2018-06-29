#################################################################################################

##              This WDL script marks the duplicates on input sorted BAMs                ##

#                              Script Options
#      -I      "Input BAM Files"                             (Required)      
#      -O      "Output BAM Files"                            (Required) 
#      -M      "File to write Duplication Metrics            (Required) 

#################################################################################################

task MarkDuplicatesTask {
   File Aligned_Sorted_Bam                  # Input Sorted BAM File
   String sampleName                               # Name of the Sample
   String Exit_Code                                # File to capture exit code
   String JAVA                                     # Variable path to Java
   String PICARD                                   # Variable path to Picard 
   String Failure_Logs                             # Variable to capture Failure Reports
   String BashScriptPath
   String dollar = "$" 
   Int Flag = 0                                    # Variable to Notify user of completion of Alignment Block
   
   command {
   
      /bin/bash ${BashScriptPath} ${Aligned_Sorted_Bam} ${sampleName} ${Exit_Code} ${JAVA} ${PICARD} ${Failure_Logs}

   }
   
   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs      
   output {
      File Aligned_Sorted_Dedupped_Bam = "${sampleName}.aligned.sorted.dedupped.bam"
      File PicardMetrics = "${sampleName}.PicardMetrics"

      #Variable to Notify user of completion of Alignment Block
      Int Notify_EndofAlignment = "${Flag}"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true

   }

}

