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
   String BashScriptPath

   command <<<
 
      /bin/bash ${BashScriptPath} ${Input_Read1} ${Input_Read2} ${sampleName} ${RefFasta} ${Ref_Amb_File} ${Ref_Dict_File} ${Ref_Ann_File} ${Ref_Bwt_File} ${Ref_Fai_File} ${Ref_Pac_File} ${Ref_Sa_File} ${BWA} > ${sampleName}.sam
 
    >>>
   
   # The output block is where the output of the program is stored.
   output {
      File Aligned_Bam = "${sampleName}.sam"
      File out = stdout()
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      # Even if the command in task has a non zero exit code continue with the other tasks
      continueOnReturnCode: true
   }

}  # End of task block
