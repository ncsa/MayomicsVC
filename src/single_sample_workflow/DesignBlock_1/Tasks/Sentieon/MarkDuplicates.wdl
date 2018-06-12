#################################################################################################

##              This WDL script marks the duplicates on input sorted BAMs                ##

#                              Script Options
#      -I      "Input BAM Files"                             (Required)      
#      -O      "Output BAM Files"                            (Required) 
#      -M      "File to write Duplication Metrics            (Required) 

#################################################################################################

task MarkDuplicatesTask {

   File InputBam                   # Input Sorted BAM File
   File InputBamIdx                # Input Sorted Bam Index File
   String sampleName               # Name of the Sample
   String Error_Logs               # File to capture exit code
   String OutDir                   # Directory for output folder
   String Sentieon                 # Variable path to Sentieon 
   String Sentieon_License         # License Server Information
   File MarkDuplicates_Script      # Bash script which is called inside the WDL script`
   String Threads                  # Specifies the number of thread required per run
   Boolean Debug_Mode_EN           # Variable to check if Debud Mode is on or not

   command {

   /bin/bash ${MarkDuplicates_Script} -L ${Sentieon_License} -b ${InputBam} -s ${sampleName} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}

   }


   output {

      File AlignedSortedDeduppedBam = "${OutDir}/${sampleName}.deduped.bam"
      File AlignedSortedDeduppedBamIndex = "${OutDir}/${sampleName}.deduped.bam.bai"
      File DedupMetrics = "${OutDir}/${sampleName}.dedup_metrics.txt"
      String sName = sampleName
      String SentieonPath = Sentieon
      String LicenseFile = Sentieon_License
      String DebugMode = Debug_Mode_EN
      String ErrLogs = Error_Logs
      String OutputDir = OutDir
      String ThreadCount = Threads

   }
}
