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
   String TempDir                       # Temporary Directory for Picard
   String OutDir                        # Directory for output folder
   String JAVA                          # Variable path to Java
   String PICARD                        # Variable path to Picard 
   File Bash_Script                     # Bash script which is called inside the WDL script
   String Threads                       # Specifies the number of thread required per run

   # Flags to be passed to the scripts

   String b="-b"
   String s="-s"
   String T="-T"
   String O="-O"
   String J="-J"
   String P="-P"
   String t="-t"
   String e="-e"
   
   command {
   
   /bin/bash ${Bash_Script} ${b} ${InputBam} ${s} ${sampleName} ${T} ${TempDir} ${O} ${OutDir} ${J} ${JAVA} ${P} ${PICARD} ${t} ${Threads} ${e} ${Error_Logs}  
   
   }
   

   output {
   
      File AlignedSortedDeduppedBam = "${OutDir}/${sampleName}.deduped.bam"
      File AlignedSortedDeduppedBamIndex = "${OutDir}/${sampleName}.deduped.bai"
      File PicardMetrics = "${OutDir}/${sampleName}.picard.metrics"

   }
}

workflow RunMarkDuplicatesTask {

   call MarkDuplicatesTask

}


