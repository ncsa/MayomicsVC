###################################################################################################
#
###              This WDL script intiates Samtools view to convert SAM Files to BAM               ##
#
###                              Script Options
###      -@      "Number of Threads"                     (Optional)
###      -o      "Output to a File"                      (Required)
#
####################################################################################################



task Samtools {

   Array[File] Aligned_Sam             # Array of input aligned SAM Files
   String sampleName                   # Name of the Sample
   String Exit_Code                    # A File that captures the list of samples that fails the Samtools step
   String SAMTOOLS                     # Variable where path the to Samtools is defined 

   command {

      # The SAMTOOLS view is used to convert Sam Files to BAM
      ${SAMTOOLS} view -@ 17 -o ${sampleName}.aligned.bam ${sep=',' Aligned_Sam}

   }

   # The output block is where the output of the program is stored. 
   # Since there are multiple output being generated, the glob function is used to capture output
   output {
      Array[File] Aligned_Bam = glob("${sampleName}.aligned.bam")
   }

   # The runtime block is used to specify Cromwell of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true
   }
}
