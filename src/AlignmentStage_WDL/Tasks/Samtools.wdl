###################################################################################################
#
###              This WDL script intiates Samtools view to convert SAM Files to BAM               ##
#
###                              Script Options
###      -@      "Number of Threads"                     (Optional)
###      -o      "Output to a File"                      (Required)
#
####################################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentStage_WDL/Workflows/TestAlignBlock.wdl" as AlignmentB

task Samtools {

   File Bam             # Array of input aligned SAM Files
   String sampleName                   # Name of the Sample
   String Exit_Code                    # A File that captures the list of samples that fails the Samtools step
   String SAMTOOLS                     # Variable where path the to Samtools is defined 

   command {

      # The SAMTOOLS view is used to convert Sam Files to BAM
      ${SAMTOOLS} sort -o ${sampleName}.bam -n ${Bam}

   }

   # The output block is where the output of the program is stored. 
   # Since there are multiple output being generated, the glob function is used to capture output
   output {
      File Obam = "${sampleName}.bam"
   }

   # The runtime block is used to specify Cromwell of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true
   }
}

workflow Call_Samtools {
   
   call AlignmentB.AlignBlock_Run
   call Samtools {
      
      input:
         Bam = AlignBlock_Run.OBama
   }
}
