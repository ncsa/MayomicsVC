#######################################################################################################

##              This WDL script performs dort on Input BAM File                ##

##                              Script Options
##      -c      "Number of Threads"                                                      (Optional)
##      -i      "Creates a BAM Index File for the final sorted output "                  (Optional)      
##      -o      "Final Output is written to a file specified"                            (Optional) 

#######################################################################################################

# The task block is used to perform Novosort on input BAM File

task Novosort {
   Array[File] Aligned_Bam                 # Input BAM File
   String sampleName                       # Name of the Sample
   String Exit_Code                        # File to capture exit code
   String SORT                             # Variable path to Novosort
   String Failure_Logs                     # Variable to capture Failure Reports

   command {

      # Check to see if input files are non-zero
      [ -s ${sep=',' Aligned_Bam} ] || echo "Input BAM File is Empty" >> ${Failure_Logs}

      # Novosort Tools is used to created sort BAM Files 
      ${SORT} -c 36 -i -o ${sampleName}.aligned.sorted.bam ${sep=',' Aligned_Bam}
      
      # The 'if' check to see if any of the samples have failed this step          
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the Novosort Step' >> ${Exit_Code}
      fi
   }

   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs   
   output {  
      Array[File] Aligned_Sorted_Bam = glob("${sampleName}.aligned.sorted.bam")
      
      # Variable to Email user at the end of step in workflow
      Int Email_Flag = 2
   }
   
   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {  
      continueOnReturnCode: true
   }
} # End of Task Block
