#########################################################################################################

###       This WDL script is performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "DesignBlock_1/Tasks/Sentieon/alignment.wdl" as ALIGN

workflow CallalignmentTask {

   # Call to the ReadMappingTask
   call ALIGN.alignmentTask

}
