#################################################################################################

task MarkDuplicatesTask {
   File Aligned_Sorted_Bam                  # Input Sorted BAM File
   String JAVA                                     # Variable path to Java
   String PICARD                                   # Variable path to Picard 
   String sampleName
   
   command {
     
      ${JAVA} -Xmx2g -jar ${PICARD} MarkDuplicates I=${Aligned_Sorted_Bam} O=${sampleName}.aligned.sorted.dedupped.bam M=${sampleName}.PicardMetrics ASSUME_SORTED=true CREATE_INDEX=true
         
   }
   
   output {
      File Aligned_Sorted_Dedupped_Bam = "${sampleName}.aligned.sorted.dedupped.bam"
      #File metrics = "${sampleName}.PicardMetrics"

   }

   runtime {
      continueOnReturnCode: true

   }

}

workflow CallMarkDuplicatesTask {

   String JAVA                               
   String PICARD                           

   Array[File] inputAlignedBam
   Array[File] GlobalsampleName

   scatter(AlignedBam in inputAlignedBam) {
      
      call MarkDuplicatesTask {
         input:
            sampleName = AlignedBam.left, 
            Aligned_Sorted_Bam = AlignedBam.right,
            JAVA = JAVA,
            PICARD = PICARD
      }
      
   }

   output {
   Array[File] Global_Aligned_Sorted_Dedupped_Bam = MarkDuplicatesTask.Aligned_Sorted_Dedupped_Bam
   #Array[File] Global_metrics = MarkDuplicatesTask.metrics
    
   }      

}
