#######################################################################################################

##              This WDL script performs sort on Input BAM File       

##                              Script Options
##      -c      "Number of Threads"                                                      (Optional)
##      -i      "Creates a BAM Index File for the final sorted output "                  (Optional)      
##      -o      "Final Output is written to a file specified"                            (Optional) 

#######################################################################################################

# The task block is used to perform Novosort on input BAM File

task NovosortTask {
   File Aligned_Bam                 # Input BAM File
   String sampleName                       # Name of the Sample
   String SORT                             # Variable path to Novosort


   command {


      # Novosort Tools is used to created sort BAM Files 
      ${SORT} -c 36 -m 20G -i -o ${sampleName}.aligned.sorted.bam ${Aligned_Bam}

}
      
   output {  
      File Aligned_Sorted_Bam = "${sampleName}.aligned.sorted.bam"
      String Global_sampleName = "${sampleName}"
   }
   
   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {  
      continueOnReturnCode: true
   }
}


workflow CallNovosortTask {

   String SORT 
   #Array[Pair[String, File]] inputAlignedBam1	
   Array[Array[File]] inputAlignedBam

   scatter(AlignedBam in inputAlignedBam) {


      call NovosortTask {
         input :
            SORT = SORT,
            sampleName = AlignedBam[0],
            Aligned_Bam = AlignedBam[1]
      }
  }

   output {
      Array[File] Global_Aligned_Sorted_Bam = NovosortTask.Aligned_Sorted_Bam
      #Array[Pair[String, File]] Global_out = inputAlignedBam
   }

}
