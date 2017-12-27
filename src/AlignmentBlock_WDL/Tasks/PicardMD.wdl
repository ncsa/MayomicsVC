#################################################################################################

##              This WDL script marks the duplicates on input sorted BAMs                ##

#                              Script Options
#      -I      "Input BAM Files"                             (Required)      
#      -O      "Output BAM Files"                            (Required) 
#      -M      "File to write Duplication Metrics            (Required) 

#################################################################################################

task Picard_MarkDuplicates {
   Array[File] Aligned_Sorted_Bam                  # Input Sorted BAM File
   String sampleName                               # Name of the Sample
   String Exit_Code                                # File to capture exit code
   String JAVA                                     # Variable path to Java
   String PICARD                                   # Variable path to Picard 
   String Failure_Logs                             # Variable to capture Failure Reports
   String dollar = "$" 
   Int Flag = 1000
   
   command {
     
      # Check to see if input files are non-zero
      [ -s ${sep=',' Aligned_Sorted_Bam} ] || echo "Aligned Sorted Bam File is Empty" >> ${Failure_Logs}
 
      # Picard Mark Duplicates is used to mark duplicates on input sorted BAMs
      ${JAVA} -Xmx2g -jar ${PICARD} MarkDuplicates I=${sep=',' Aligned_Sorted_Bam } O=${sampleName}.aligned.sorted.dedupped.bam M=${sampleName}.PicardMetrics ASSUME_SORTED=true CREATE_INDEX=true
         
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the Mark Duplicates Step' >> ${Exit_Code}
      fi

      [ ! -f ${sampleName}.aligned.sorted.dedupped.bam ] && echo "aligned sorted dedupped bam not created" >> ${Failure_Logs}
   }
   
   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs      
   output {
      Array[File] Aligned_Sorted_Dedupped_Bam = glob("${sampleName}.aligned.sorted.dedupped.bam")
      Array[File] PicardMetrics = glob("${sampleName}.PicardMetrics")

      #Variable to Notify user of completion of Alignment Block
      Int Notify_EndofAlignment = "${Flag}"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true

   }

}

