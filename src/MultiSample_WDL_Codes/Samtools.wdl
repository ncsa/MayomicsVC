###################################################################################################

##		This WDL script intiates Samtools view to convert SAM Files to BAM		 ##

##                              Script Options
##	-@	"Number of Threads" 			(Optional)
##	-o	"Output to a File"			(Required)

###################################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/MultiSample_WDL_Codes/BWA_Mem.wdl" as BWA




# The Task block is where the variables and the functions are defined for performing a certain task
task Samtools {

   Array[File] Aligned_Sam		# Array of input aligned SAM Files
#   String sampleName			# Name of the Sample
#   String Exit_Code			# A File that captures the list of samples that fails the Samtools step
   String SAMTOOLS			# Variable where path the to Samtools is defined 

   command {
  
      # The SAMTOOLS tool is used to convert Sam Files to BAM
      SAMTOOLS view -@ 17 -o aligned.bam ${sep=',' Aligned_Sam}

   }

   # The output block is where the output of the program is stored. 
   # Since there are multiple output being generated, the glob function is used to capture output
   output {
     Array[File] Aligned_Bam = glob("aligned.bam")
   }

   # The runtime block is used to specify Cromwell of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true
   }
}

workflow Samtools_Run {
 
   call BWA.BWA_Mem_Run

   scatter(sample in BWA_Mem_Run.inputsamples) {
      call Samtools {
         input :
            Aligned_Sam = BWA_Mem_Run.Aligned_Sam,
      }
   }  
}



