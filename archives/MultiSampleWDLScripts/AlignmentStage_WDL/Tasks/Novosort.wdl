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


   command {

      # The 'set' command is used set and examine shell options, as well as set positional parameters
      set -x

      # Check to see if input files are non-zero
      [ -s ${Aligned_Bam} ] || echo "Input BAM File is Empty" >> ${Failure_Logs}

      # Record the start of the program execution
      StartTime=`date +%s`

      # Novosort Tools is used to created sort BAM Files 
      ${NOVOSORT} -c 36 -i -o ${sampleName}.aligned.sorted.bam ${sep=',' Aligned_Bam}

      # Record the End of the program execution
      EndTime=`date +%s`

      # Calculate the time to complete the run BWA and Samtools on each sample
      echo "${sampleName} ran Novosort for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
      
      # The 'if' check to see if any of the samples have failed this step          
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the Novosort Step' >> ${Exit_Code}
      fi

      [ ! -f ${sampleName}.aligned.sorted.bam ] && echo "aligned sorted bam not created" >> ${Failure_Logs}
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
