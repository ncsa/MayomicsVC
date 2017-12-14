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

task BWA_Mem {
        
   File Input_Read1		# Input Read File		 (REQUIRED)
   File Input_Read2		# Input Read File		 (Optional)
   String sampleName		# Name of the Sample
   File RefFasta		# Reference FASTA file		 (REQUIRED)

   File Ref_Amb_File		#
   File Ref_Dict_File		#
   File Ref_Ann_File		#
   File Ref_bwt_File		# These are reference files that are provided as implicit inputs
   File Ref_fai_File		# to the WDL Tool to help perform the alignment
   File Ref_pac_File		#
   File Ref_sa_File		#

   String BWA			# Variable path to BWA MEM Tool
   String SAMTOOL               # variable path to Samtools
   String Exit_Code		# Variable capture exit code
   String Failure_Logs          # Variable to capture Failure messages

   command <<<

      # pipefall command sets the exit status to the exit code of the last program to exit non-zero 
      # or zero if all exited successfully    
      set -o pipefail
      
      # To check if the executables are present
      
      [ ! -f ${BWA} ] && echo "BWA does not exist" >>${Failure_Logs}

      [ ! -f ${SAMTOOL} ] && echo "Samtools does not exist" >>${Failure_Logs}

      [ ! -f ${RefFasta} ] && echo "Reference File does does not exist" >>${Failure_Logs}

      # To check if the Input Files are present

      [ -s ${RefFasta} ] || echo "Reference File is Empty" >>${Failure_Logs}

      [ -s ${Input_Read1} ] || echo "Input Read 1 is Empty" >>${Failure_Logs}

      [ -s ${Input_Read2} ] || echo "Input Read 2 is Empty" >>${Failure_Logs}

      # BWA Mem Tool is used to create aligned SAM file from the input FASTA File
      ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} | ${SAMTOOL} view -@ 2 -bSu -> ${sampleName}.aligned.bam
      
      
      # The 'if' check to see if any of the samples have failed this step
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the BWA Mem Step ' >> ${Exit_Code}
      fi
  
    >>>
   
   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs
   output {
      Array[File] Aligned_Bam = glob("${sampleName}.aligned.bam")
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      # Even if the command in task has a non zero exit code continue with the other tasks
      continueOnReturnCode: true
      failOnStderr: false
   }

}  # End of task block
