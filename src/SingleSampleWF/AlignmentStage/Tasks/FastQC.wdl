##########################################################################################################

##            This WDL script performs Quality Control on input FastQ files            ##

##                                      Script Options

##      -t             "Number of Threads"                                                   (Optional)
##      -j             "Full path to the java binary used to launch fastqc"                  (Optional)      
##      -o             "All output files created in the specified output directory"          (Optional) 
##      --extract      "Zipped output file uncompressed in the same directory"               (Optional) 

##########################################################################################################

task FastqQualityControlTask {

   File Input_Read1             # Input Read File                (REQUIRED)
   File Input_Read2             # Input Read File                (Optional)
   String sampleName            # Name of the Sample
   String FastQCDir             # Output Folder directory where the HTML report and the summary text is stored
  
   String FASTQC                # Variable path for FastQC
   String JAVA                  # Variable path for Java
   String Failure_Logs          # Variabkle to capture Failure Reports
   String Exit_Code             # Variable capture exit code

   command <<<

      # Check to see if input files are non-zero
      [ -s ${Input_Read1} ] || echo "Input Read1 FastQ is Empty" >> ${Failure_Logs}
      [ -s ${Input_Read2} ] || echo "Input Read2 FastQ is Empty" >> ${Failure_Logs} 

      #FastQC takes a FastQ file and runs a series of tests on it to generate a comprehensive QC report
      ${FASTQC} --extract -j ${JAVA} -o ${FastQCDir} ${Input_Read1} ${Input_Read2}

      if [ $? -ne 0 ]; then
         echo "${sampleName} has failed at the FASTQ/BAM File Quality Control Step" >> ${Exit_Code}
      fi

      [ ! -d ${FastQCDir} ] && echo "FASTQC directory has not been created" >> ${Failure_Logs}
   
   >>>

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for calls
   runtime {
      # Even if the command in task has a non zero exit code continue with the other tasks
      continueOnReturnCode: true
   }

}  # End of task block
     
