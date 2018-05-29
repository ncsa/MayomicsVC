#########################################################################################################

###       This WDL script is performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "DesignBlock_1/Tasks/BWAMemSamtoolsView.wdl" as BWAMEMSAMTOOLSVIEW

workflow CallReadMappingTask {

   # Call to the ReadMappingTask
   call BWAMEMSAMTOOLSVIEW.ReadMappingTask

}
