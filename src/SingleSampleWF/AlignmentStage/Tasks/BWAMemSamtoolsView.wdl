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
   File REF_GENOME		   # Reference FASTA file        (REQUIRED)
       
   File Ref_Amb_File		   #
   File Ref_Dict_File		   #
   File Ref_Ann_File		   #
   File Ref_Bwt_File		   # These are reference files that are provided as implicit inputs
   File Ref_Fai_File		   # to the WDL Tool to help perform the alignment
   File Ref_Pac_File		   #
   File Ref_Sa_File		   #
        
   String BWA			   # Variable path to BWA MEM Tool
   String SAMTOOLS                  # variable path to Samtools
   String Exit_Code		   # Variable capture exit code
   String Failure_Logs             # Variable to capture Failure reports
   String dollar = "$"
        
   command <<<
        
      # The 'set' command is used to set and examine shell options, as well as set positional parameters
      set -x
       
      # Check to see if input files are non-zero
      [[ -s ${Input_Read1} ]] || echo -e "Input Read 1 File is Empty" >> ${Failure_Logs} 
      [[ -s ${Input_Read2} ]] || echo -e "Input Read 2 File is Empty" >> ${Failure_Logs}
       
      # Record the start of the program execution
      StartTime=`date +%s`
          
      # BWA Mem Tool is used to create aligned SAM file from the input FASTA File
      # PIPESTATUS is an internal bash variable which holds the exit code of commands in the pipe
      ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${REF_GENOME} ${Input_Read1} ${Input_Read2} | ${SAMTOOLS} view -@ 17 -bSu -> ${sampleName}.aligned.bam; B=(${dollar}{PIPESTATUS[*]})
          
      # Record the End of the program execution
      EndTime=`date +%s`
          
      # Calculate the time to complete the run BWA and Samtools on each sample
      echo -e "${sampleName} ran BWA Mem and Samtools View for ${dollar}((${dollar}{EndTime} - ${dollar}{StartTime})) seconds" >> ${Failure_Logs}
           
      # Checks to see if either BWA or Samtools has a non-zero exit code
      if [ ${dollar}{B[0]} -ne 0 ] 
      then
         echo -e "${sampleName} exited BWA with code ${dollar}{B[0]}" >> ${Exit_Code}
      fi
        
      if [ ${dollar}{B[1]} -ne 0 ]
      then
         echo -e "${sampleName} exited SAMTOOLS with code ${dollar}{B[1]}" >> ${Exit_Code}
      fi
          
      # Check to see if the BAM file is created or not 
      [[ ! -f ${sampleName}.aligned.bam ]] && echo -e "aligned bam not created" >> ${Failure_Logs}
         
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
