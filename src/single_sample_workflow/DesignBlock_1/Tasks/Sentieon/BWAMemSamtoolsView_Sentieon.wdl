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

   String Error_Logs               # File Path to ErrorLogs
   String OutDir                   # Directory for output folder
   String Threads                  # Specifies the number of thread required per run
   File Bash_Script                # Bash script which is called inside the WDL script

   # Flags to be passed to the scripts
 
   String SE="-SE"
   String r="-r"
   String R="-R"
   String s="-s"
   String G="-G"
   String O="-O"
   String t="-t"
   String e="-e"
   String g="-g"
   String p="-p"

   command {

      # Check to see if the Input FastQ is Singled Ended or not
      if [[ ${Is_Single_End} == false ]] 
      then
         /bin/bash ${Bash_Script} ${SE} ${Is_Single_End} ${g} ${Group} ${r} ${Input_Read1} ${R} ${Input_Read2} ${s} ${sampleName} ${p} ${Platform} ${G} ${RefFasta} ${O} ${OutDir} ${S} ${Sentieon} ${t} ${Threads} ${e} ${Error_Logs}

      else
         /bin/bash ${Bash_Script} ${SE} ${Is_Single_End} ${g} ${Group} ${r} ${Input_Read1} ${s} ${sampleName} ${p} ${Platform} ${G} ${RefFasta} ${O} ${OutDir} ${S} ${Sentieon} ${t} ${Threads} ${e} ${Error_Logs}
      fi

   }

   # The output block is where the output of the program is stored
   output {

      File OutSam = "${OutDir}/${sampleName}.sam"
      File OutBam = "${OutDir}/${sampleName}.bam"
      File SortBam = "${OutDir}/${sampleName}.sorted.bam"

   }

} 


workflow CallReadMappingTask {

call ReadMappingTask

}
