#################################################################################################

##              This WDL script marks the duplicates on input sorted BAMs                ##

#                              Script Options
#      -I      "Input BAM Files"                             (Required)      
#      -O      "Output BAM Files"                            (Required) 
#      -M      "File to write Duplication Metrics            (Required) 

#################################################################################################

task MarkDuplicatesTask {

   File InputBam                        # Input Sorted BAM File
   String sampleName                    # Name of the Sample
   String Error_Logs                    # File to capture exit code
   String OutDir                        # Directory for output folder
   String Sentieon                      # Variable path to Sentieon 
   File Bash_Script                     # Bash script which is called inside the WDL script`
   String Threads                       # Specifies the number of thread required per run

   command {

   /bin/bash ${Bash_Script} -b ${InputBam} -s ${sampleName} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs}

   }


   output {

      File AlignedSortedDeduppedBam = "${OutDir}/${sampleName}.deduped.bam"
      File AlignedSortedDeduppedBamIndex = "${OutDir}/${sampleName}.deduped.bam.bai"
      File DedupMetrics = "${OutDir}/${sampleName}.dedup_metrics.txt"

   }
}
