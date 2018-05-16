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
   String BashScriptPath

   command <<<
   
   /bin/bash ${BashScriptPath} ${Input_Read1} {Input_Read2} ${sampleName} ${FastQCDir} ${FASTQC} $JAVA{} ${Failure_Logs} ${Exit_Code}
   
   >>>

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for calls
   runtime {
      # Even if the command in task has a non zero exit code continue with the other tasks
      continueOnReturnCode: true
   }

}  # End of task block
     
