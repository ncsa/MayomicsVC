#######################################################################################################

##              This WDL script performs sort on Input BAM File       

##                              Script Options
##      -c      "Number of Threads"                                                      (Optional)
##      -i      "Creates a BAM Index File for the final sorted output "                  (Optional)      
##      -o      "Final Output is written to a file specified"                            (Optional) 

#######################################################################################################

# The task block is used to perform Novosort on input BAM File

task NovosortTask {
   File Aligned_Bam                 # Input BAM File
   String sampleName                       # Name of the Sample
   String Exit_Code                        # File to capture exit code
   String NOVOSORT                             # Variable path to Novosort
   String Failure_Logs                     # Variable to capture Failure Reports
   String BashScriptPath


   command {

      /bin/bash ${BashScriptPath} ${Aligned_Bam} ${sampleName} ${Exit_Code} ${Novosort} ${Failure_Logs}

   }

   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs   
   output {  
      File Aligned_Sorted_Bam = "${sampleName}.aligned.sorted.bam"
   }
   
   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {  
      continueOnReturnCode: true
   }
} # End of Task Block
