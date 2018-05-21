#########################################################################################################

###       This WDL script is performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "AlignmentStage_WDL/Tasks/BWAMemSamtoolView.wdl" as BWAMEMSAMTOOLVIEW

workflow CallReadMappingTask {

      call BWAMEMSAMTOOLVIEW.ReadMappingTask

} # End of Workflow block

