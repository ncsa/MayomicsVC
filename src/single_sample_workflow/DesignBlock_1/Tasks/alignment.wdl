###########################################################################################

##              This WDL script performs alignment using BWA Mem                ##

##                              Script Options
#       -t        "Number of Threads"                         (Optional)
#       -SE       "Single Ended Reads specification"          (Required)
#       -r        "Left Fasta File"                           (Required)
#       -R        "Right Fasta File"                          (Optional)
#       -s        "Name of the sample"                        (Optional)
#       -S        "Path to the Sentieon Tool"                 (Required)
#       -O        "Directory for the Output"                  (Required)
#       -L        "Sentieon License File"                     (Required)
#       -e        "Path to the Error Log File"                (Required)
#       -d        "Debug Mode Toggle"                         (Optional)
#       -g        "Group"                                     (Required)
#       -p        "Platform"                                  (Required)

###########################################################################################

# The Task block is where the variables and the functions are defined for performing a certain task

task alignmentTask {

   File RefFasta                   # Reference Input Fasta File
   File Input_Read1                # Input Read File             (REQUIRED)
   String Input_Read2              # Input Read File             (Optional)
   String sampleName               # Name of the Sample

   File Ref_Amb_File               #
   File Ref_Dict_File              #
   File Ref_Ann_File               #
   File Ref_Bwt_File               # These are reference files that are provided as implicit inputs
   File Ref_Fai_File               # to the WDL Tool to help perform the alignment
   File Ref_Pac_File               #
   File Ref_Sa_File                #

   String Sentieon_License         # Sentieon License server
   String Sentieon                 # Path to Sentieon
   String Group
   String Platform

   Boolean Is_Single_End           # Variable to check if single ended or not
   Boolean Debug_Mode_EN           # Variable to check if Debud Mode is on or not
   String Error_Logs               # File Path to ErrorLogs
   String OutDir                   # Directory for output folder
   String Threads                  # Specifies the number of thread required per run
   File Alignment_Script           # Bash script which is called inside the WDL script

   command {

      # Check to see if the Input FastQ is Singled Ended or not
      if [[ ${Is_Single_End} == false ]] 
      then
         /bin/bash ${Alignment_Script} -L ${Sentieon_License} -P ${Is_Single_End} -g ${Group} -r ${Input_Read1} -R ${Input_Read2} -s ${sampleName} -p ${Platform} -G ${RefFasta} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}

      else
         /bin/bash ${Alignment_Script} -L ${Sentieon_License} -P ${Is_Single_End} -g ${Group} -r ${Input_Read1} -R ${sampleName} -P ${Platform} -G ${RefFasta} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}
      fi

   }

   # The output block is where the output of the program is stored
   output {

      File OutSam = "${OutDir}/${sampleName}.sam"
      File SortBam = "${OutDir}/${sampleName}.sorted.bam"
      File SortBamIdx = "${OutDir}/${sampleName}.sorted.bam.bai"
      String sName= sampleName
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

   }

} 

