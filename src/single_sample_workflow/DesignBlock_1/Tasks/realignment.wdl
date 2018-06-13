###########################################################################################

##              This WDL script performs realignment using Sentieon                ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -b        "Input Deduped Bam"                                     (Required)
#       -k        "List of Known Sites"                                   (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -O        "Directory for the Output"                              (Required)
#       -L        "Sentieon License File"                                 (Required)
#       -e        "Path to the Error Log File"                            (Required)
#       -d        "Debug Mode Toggle"                                     (Optional)

###########################################################################################

task realignmentTask {

   File InputBam                   # Input Sorted Deduped Bam
   File InputBamIdx                # Input Sorted Deduped Bam Index
   File RefFasta                   # Reference Genome
                                   
   File Ref_Amb_File               #
   File Ref_Dict_File              #
   File Ref_Ann_File               #
   File Ref_Bwt_File               # These are reference files that are provided as implicit inputs
   File Ref_Fai_File               # to the WDL Tool to help perform the realignment
   File Ref_Pac_File               #
   File Ref_Sa_File                #
                                   
   String sampleName               # Name of the Sample

   File Known_Sites                # List of known sites
   File KnownSitesIdx              # Index file for the known sites
   String Threads                  # No of Threads for the Tool
   String Sentieon_License         # Sentieon License Information
   String Sentieon                 # Path to Sentieon
   String OutDir                   # Output Directory
   Boolean Debug_Mode_EN           # Enable or Disable Debug Mode
   String Error_Logs               # Filepath to Error Logs
   
   File Realignment_Script         # Path to bash script called within WDL script
 
   command {

      /bin/bash ${Realignment_Script} -L ${Sentieon_License} -s ${sampleName} -b ${InputBam} -G ${RefFasta} -k ${Known_Sites} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}

   }

  
   # The output block is where the output of the program is stored
   output {
  
      File OutBam = "${OutDir}/${sampleName}.realigned.bam"
      File OutBamIdx = "${OutDir}/${sampleName}.realigned.bam.bai"
      String SentieonPath = Sentieon
      String LicenseFile = Sentieon_License
      Boolean DebugMode = Debug_Mode_EN
      String ErrLogs = Error_Logs
      String OutputDir = OutDir
      String ThreadCount = Threads

      File FastaRef = RefFasta
      File RefAmbFile = Ref_Amb_File
      File RefDictFile = Ref_Dict_File
      File RefAnnFile = Ref_Ann_File
      File RefBwtFile = Ref_Bwt_File
      File RefFaiFile = Ref_Fai_File
      File RefPacFile = Ref_Pac_File
      File RefSaFile = Ref_Sa_File
      File KnownSites = Known_Sites
      File KnownSitesIDX = KnownSitesIdx
       
   }  
   
} 
