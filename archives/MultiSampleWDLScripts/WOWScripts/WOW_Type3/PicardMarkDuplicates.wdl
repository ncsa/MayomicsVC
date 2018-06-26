#################
import /projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/WOWScripts/WOW_Type3/MasterWF.wdl

task MarkDuplicatesTask {
   File Aligned_Sorted_Bam                  # Input Sorted BAM File
   String sampleName                               # Name of the Sample
   String JAVA                                     # Variable path to Java
   String PICARD                                   # Variable path to Picard 
   
   command {
     
 
      # Picard Mark Duplicates is used to mark duplicates on input sorted BAMs
      ${JAVA} -Xmx2g -jar ${PICARD} MarkDuplicates I=${Aligned_Sorted_Bam} O=${sampleName}.aligned.sorted.dedupped.bam M=${sampleName}.PicardMetrics ASSUME_SORTED=true CREATE_INDEX=true
         
   }
   
   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs      
   output {
      File Aligned_Sorted_Dedupped_Bam = "${sampleName}.aligned.sorted.dedupped.bam"
      File PicardMetrics = "${sampleName}.PicardMetrics"

      #Variable to Notify user of completion of Alignment Block
      Int Notify_EndofAlignment = "${Flag}"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true

   }

}

workflow CallMarkDuplicatesTask {
   
      
   
   scatter(bam in AlignedSortedBam) {

      call MarkDuplicatesTask {
         input :
            JAVA = JAVA,
            PICARD = PICARD,
            sampleName = bam[0],
            Aligned_Sorted_Bam = bam[1]

      }
   }
 
 
}
