###########################

#Trim Input Sequences using Cutadapt

task TrimInputSequencesTask {

   File Input_Read1                # Input Read File             (REQUIRED)
   String Input_Read2              # Input Read File             (Optional)
   File Adapters                   # Adapter FastA File          (REQUIRED) 
   String OutDir                   # Directory for output folder
   String CutAdapt                 # Path to CutAdapt Tool
   String Threads                  # Specifies the number of thread required per run
   Boolean Is_Single_End           # Variable to check if single ended or not
   String Error_Logs               # File Path to ErrorLogs
   File Bash_Script                # Bash script which is called inside the WDL script
   String sampleName               # Name of the Sample

   command {

      # Check to see if the Input FastQ is Singled Ended or not
      if [[ ${Is_Single_End} == false ]]
      then
         /bin/bash ${Bash_Script} -SE ${Is_Single_End} -r ${Input_Read1} -R ${Input_Read2} -s ${sampleName} -A ${Adapters} -O ${OutDir} -C ${CutAdapt} -t ${Threads} -e ${Error_Logs}

      else
         /bin/bash ${Bash_Script} -SE ${Is_Single_End} -r ${Input_Read1} -s ${sampleName} -A ${Adapters} -O ${OutDir} -C ${CutAdapt} -t ${Threads} -e ${Error_Logs}
      fi

   }

   # The output block is where the output of the program is stored
   output {

      File TrimmedInputRead = "${Outdir}/${sampleName}.trimmed.fq.gz"
      File TrimmedInputRead1 = "${Outdir}/${sampleName}.read1.trimmed.fq.gz"
      File TrimmedInputRead2 = "${Outdir}/${sampleName}.read2.trimmed.fq.gz"

   }

}

