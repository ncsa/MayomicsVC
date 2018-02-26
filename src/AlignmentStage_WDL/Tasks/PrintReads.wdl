###########################################################################################
#
##            This WDL script writes out sequencing data using Print Reads               ##

##                                  Script Options                                     ##

##             -T                "Name of the tool to run"                  (Required)
###             --out            "Write output to this BAM filename"        (Required)
###             -nct             "Number of Threads"                        (Optional)
###             -R               "Reference FASTA file"                     (Optional)
############################################################################################

task PrintReads {
   File RefFasta
   File Aligned_Sorted_Dedupped_Bam
   File Recal_Report
   String sampleName
   String JAVA
   String GATK
   String Exit_Code
   String Failure_Logs
   String dollar = "$"         

   command {
      
      # Check to see if input files are non-zero
      [ -s ${Aligned_Sorted_Dedupped_Bam} ] || echo "Aligned Sorted Dedupped Bam File is Empty" >> ${Failure_Logs}


      ${JAVA} -Xmx16g -jar ${GATK} -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Bam} -T PrintReads -BQSR ${Recal_Report} --out ${sampleName}.aligned.sorted.dedupped.realigned.recalibrated.bam -nct 6
    
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the Print Reads Step' >> ${Exit_Code}
      fi

      [ ! -f ${sampleName}.aligned.sorted.dedupped.realigned.recalibrated.bam ] && echo "Aligned sorted dedupped realigned recalibrated bam not created" >> ${Failure_Logs}
   }

   output {
      File Aligned_Sorted_Dedupped_Realigned_Recalibrated_bam = "${sampeName}.aligned.sorted.dedupped.realigned.recalibrated.bam"

   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true

   }

}
