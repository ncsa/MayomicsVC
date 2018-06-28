#########################################################################################################

###       This WDL script performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "wdl_scripts/DesignBlock_1/Tasks/alignment.wdl" as ALIGN

workflow CallalignmentTask {

   # Call to the ReadMappingTask
   call ALIGN.alignmentTask

}
