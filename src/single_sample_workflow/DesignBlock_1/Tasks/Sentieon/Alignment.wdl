###########################################################################################

##              This WDL script performs alignment using BWA Mem                ##

##                              Script Options
##      -t      "Number of Threads"                                     (Optional)
##      -M      "Mark shorter split hits as secondary"                  (Optional)      
##      -k      "Minimun Seed Length"                                   (Optional) 
##      -I      "The input is in the Illumina 1.3+ read format"         (Optional) 
##      -R      "Complete read group header line"                       (Optional) 

###########################################################################################

# The Task block is where the variables and the functions are defined for performing a certain task

task ReadMappingTask {

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
         /bin/bash ${Alignment_Script} -SE ${Is_Single_End} -G ${Group} -r ${Input_Read1} -R ${Input_Read2} -s ${sampleName} -P ${Platform} -G ${RefFasta} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}

      else
         /bin/bash ${Alignment_Script} -SE ${Is_Single_End} -G ${Group} -r ${Input_Read1} -R ${sampleName} -P ${Platform} -G ${RefFasta} -O ${OutDir} -S ${Sentieon} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}
      fi

   }

   # The output block is where the output of the program is stored
   output {

      File OutSam = "${OutDir}/${sampleName}.sam"
      File OutBam = "${OutDir}/${sampleName}.bam"
      File SortBam = "${OutDir}/${sampleName}.sorted.bam"

   }

} 

