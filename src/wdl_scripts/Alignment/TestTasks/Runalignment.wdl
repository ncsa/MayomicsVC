#########################################################################################################

###       This WDL script performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "src/wdl_scripts/Alignment/Tasks/alignment.wdl" as ALIGN

workflow RunAlignmentTask {

   Array[Array[File]] InputReads   # One lane per subarray with one or two input reads
   Array[String] PlatformUnit      # One platform unit per alignment task
   Boolean PairedEnd               # Variable to check if single ended or not

   Array[Int] Indexes = range(length(InputReads))

   scatter (idx in Indexes) {

      if(PairedEnd) {
         call ALIGN.alignmentTask as ALIGN_paired {
            input:
               InputRead1=InputReads[idx][0],
               InputRead2=InputReads[idx][1],
               PlatformUnit=PlatformUnit[idx]
         }
      }

      if(!PairedEnd) {
         call ALIGN.alignmentTask as ALIGN_single {
            input:
               InputRead1=InputReads[idx][0],
               InputRead2="null",
               PlatformUnit=PlatformUnit[idx]
         }
      }
   }

   output {
      # Unify outputs from scatter and filter out null entries 
      Array[File] OutputBams = select_all(flatten([ALIGN_paired.OutputBams,ALIGN_single.OutputBams]))
      Array[File] OutputBais = select_all(flatten([ALIGN_paired.OutputBais,ALIGN_single.OutputBais]))
   }

}
