###########################################################################################

##		This WDL script performs alignment using BWA Mem		##

##				Script Options
##	-t	"Number of Threads"					(Optional)
##	-M	"Mark shorter split hits as secondary"			(Optional)	
##	-k	"Minimun Seed Length"					(Optional) 
##	-I	"The input is in the Illumina 1.3+ read format"		(Optional) 
##	-R 	"Complete read group header line"			(Optional) 

###########################################################################################

# The Task block is where the variables and the functions are defined for performing a certain task

task ReadMappingTask {
        
   File Input_Read1		   # Input Read File		 (REQUIRED)
   File Input_Read2		   # Input Read File		 (Optional)
   String sampleName		   # Name of the Sample
   File RefFasta		   # Reference FASTA file		 (REQUIRED)

   File Ref_Amb_File		   #
   File Ref_Dict_File		   #
   File Ref_Ann_File		   #
   File Ref_Bwt_File		   # These are reference files that are provided as implicit inputs
   File Ref_Fai_File		   # to the WDL Tool to help perform the alignment
   File Ref_Pac_File		   #
   File Ref_Sa_File		   #

   String BWA			   # Variable path to BWA MEM Tool
   String SAMTOOL                  # variable path to Samtools
   String Exit_Code		   # Variable capture exit code
   String Failure_Logs             # Variable to capture Failure reports
   String dollar = "$"             # Variable to access internal bash variables
   String DummyVar                 # Dummy Variable created to force sequential execution
   String Email_ID                 # Variable to hold the email address
   
   command <<<

      # Check to see if input files are non-zero
      [ -s ${Input_Read1} ] || echo "Input Read 1 File is Empty" >> ${Failure_Logs} 
      [ -s ${Input_Read2} ] || echo "Input Read 2 File is Empty" >> ${Failure_Logs}
      
      # BWA Mem Tool is used to create aligned SAM file from the input FASTA File
      # PIPESTATUS is an internal bash variable which holds the exit code of commands in the pipe
      ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} | ${SAMTOOL} view -@ 17 -bSu -> ${sampleName}.aligned.bam; B=(${dollar}{PIPESTATUS[*]})

      if [ ${dollar}{B[0]} -ne 0 ] 
      then
         echo "${sampleName} exited BWA with code ${dollar}{B[0]}" >> ${Exit_Code}
         echo "${sampleName} exited BWA with code ${dollar}{B[0]}" | mailx -s "Sample Failed BWA Step" ${Email_ID}
      else
         echo "NO ERRORS" >> ${Failure_Logs}
      fi

      if [ ${dollar}{B[1]} -ne 0 ]
      then
         echo "${sampleName} exited SAMTOOLS with code ${dollar}{B[1]}" >> ${Exit_Code}
         echo "${sampleName} exited BWA with code ${dollar}{B[1]}" | mailx -s "Sample Failed SAMTOOLS Step" ${Email_ID}
      else
         echo "NO ERRORS" >> ${Failure_Logs}
      fi
     
      #The 'if' check to see if any of the samples have failed this step
      ###if [ $? -ne 0 ]; then
      ###echo '${sampleName} has failed at the BWA Mem Step ' >> ${Exit_Code}
      ###fi
     
      [ ! -f ${sampleName}.aligned.bam ] && echo "aligned bam not created" >> ${Failure_Logs}
  
    >>>
   
   # The output block is where the output of the program is stored.
   output {
      File Aligned_Bam = "${sampleName}.aligned.bam"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      # Even if the command in task has a non zero exit code continue with the other tasks
      continueOnReturnCode: true
   }

}  # End of task block
