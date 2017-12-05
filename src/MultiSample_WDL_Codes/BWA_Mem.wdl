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

   String BWA			# Variable where path the to BWA MEM Tool is defined
   String Exit_Code		# Variable used to capture the exit code


   command {

      # BWA Mem Tool is used to create aligned SAM file from the input FASTA File
      ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} > ${sampleName}.aligned.sam

      # The 'if' block copies the list of all the samples that fails the BWA Mem Step
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the BWA Mem Step ' >> ${Exit_Code}
      fi
   }
   
   # The output block is where the output of the program is stored. Since there are multiple output being generated, we use to glob function to capture the output of various samples in seperate files
   output {
      Array[File] Aligned_Sam = glob("${sampleName}.aligned.sam")
   }

   # The runtime block is used to specify the Cromwell Execution Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true
   }

}  # End of task block




# Inside the Workflow block is where the Task block is called. 
workflow BWA_Mem_Run {
   # The variable below represents a text file which has sample input information
   # Each line represents the various inputs for an individual sample
   # These inputs for a single sample are separated from one another by TAB spacing
	
   File InputSamplesFile 
   File Ref_Amb_File            #
   File Ref_Dict_File           #
   File Ref_Ann_File            #
   File Ref_bwt_File            # These are reference files that are provided as implicit inputs
   File Ref_fai_File            # to the WDL Tool to help perform the alignment
   File Ref_pac_File            #
   File Ref_sa_File 
   File RefFasta
   String BWA      
   String Exit_Code
   # The variable below reads a TSV file and store the data in an array	

   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   # The Call block is where the task to be performed is called
   
   

   # Since BWA Mem is performed on multiple samples we have to perform BWA Mem  on each sample
   # The Scatter function helps parallelize a series of identical tasks but give them different inputs
   scatter(sample in inputsamples) {
	
      call BWA_Mem {
         # The Input section inside the Call block is where input variables are assigned either with values from an array or to an output of another task 
         # This value is input for the task being called by the call function    
         input :
            # Samples 0, 1 & 2 represent information for an individual sample
            sampleName = sample[0],         
            Input_Read1 = sample[1],         
            Input_Read2 = sample[2],

            Ref_Amb_File = Ref_Amb_File,
            Ref_Dict_File = Ref_Dict_File,
            Ref_Ann_File =  Ref_Ann_File,
            Ref_bwt_File = Ref_bwt_File,
            Ref_fai_File = Ref_fai_File,
            Ref_pac_File = Ref_pac_File,
            Ref_sa_File = Ref_sa_File,
            RefFasta = RefFasta,
            BWA = BWA,
            Exit_Code = Exit_Code
      }
               
   }

   output {
      Array[Array[File]] Aligned_Sam = BWA_Mem.Aligned_Sam
      Array[String] sampleName = BWA_Mem.sampleName
   }

}  # End of Workflow block
