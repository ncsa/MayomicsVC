###############################################################################################

#       This WDL script is performs PICARD Mark Duplicates on input sorted BAM Files     ##

###############################################################################################

import "AlignmentStage_WDL/Tasks/PicardMarkDuplicates.wdl" as PICARDMARKDUPLICATES

workflow CallMarkDuplicatesTask {

      call PICARDMARKDUPLICATES.MarkDuplicatesTask 
}
