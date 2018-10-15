#########################################################################################################

###       This WDL script performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "src/wdl_scripts/Alignment/Tasks/alignment.wdl" as ALIGN

workflow CallalignmentTask {

   # Call to the ReadMappingTask
   call ALIGN.alignmentTask

}
