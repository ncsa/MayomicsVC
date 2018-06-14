###########################################################################################

##              This WDL scripts trim the Inputs Fasta File using CutAdapt               ##

##                                    Script Options                     
#         -t        "Number of Threads"                         (Required)
#         -P        "Single Ended Reads specification"          (Required)
#         -r        "Left Fastq File"                           (Required)
#         -R        "Right Fastq File"                          (Optional)
#         -s        "Name of the sample"                        (Optional)
#         -A        "Adapter File for CutAdapt"                 (Required)
#         -O        "Directory for the Output"                  (Required)
#         -C        "Path to CutAdapt Tool"                     (Required)
#         -e        "Path to the Error Log File"                (Required)
#         -d        "Debug Mode Toggle"                         (Optional)

###########################################################################################         

task trimsequencesTask {

   File Input_Read1                # Input Read File             (REQUIRED)
   String Input_Read2              # Input Read File             (Optional)
   File Adapters                   # Adapter FastA File          (REQUIRED) 
   String OutDir                   # Directory for output folder
   String CutAdapt                 # Path to CutAdapt Tool
   String Threads                  # Specifies the number of thread required per run
   Boolean Is_Single_End           # Variable to check if single ended or not
   Boolean Debug_Mode_EN           # Variable to check if Debud Mode is on or not
   String Error_Logs               # File Path to ErrorLogs
   File CutAdaptTrimming_Script    # Bash script which is called inside the WDL script
   String sampleName               # Name of the Sample

   command {

      # Check to see if the Input FastQ is Singled Ended or not
      if [[ ${Is_Single_End} == false ]]
      then
         /bin/bash ${CutAdaptTrimming_Script} -P ${Is_Single_End} -l ${Input_Read1} -r ${Input_Read2} -s ${sampleName} -A ${Adapters} -O ${OutDir} -C ${CutAdapt} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}

      else
         /bin/bash ${CutAdaptTrimming_Script} -P ${Is_Single_End} -l ${Input_Read1} -s ${sampleName} -A ${Adapters} -O ${OutDir} -C ${CutAdapt} -t ${Threads} -e ${Error_Logs} -d ${Debug_Mode_EN}
      fi

   }

   # The output block is where the output of the program is stored
   output {

      File TrimmedInputRead1 = "${OutDir}/${sampleName}.read1.trimmed.fq.gz"
      File TrimmedInputRead2 = "${OutDir}/${sampleName}.read2.trimmed.fq.gz"
      Boolean Debug_Mode = Debug_Mode_EN
      String sName = sampleName
      String ErrLogs = Error_Logs
      Boolean SingleEnd = Is_Single_End
      String ThreadCount = Threads
      String OutputDir = OutDir
      
   }

}

