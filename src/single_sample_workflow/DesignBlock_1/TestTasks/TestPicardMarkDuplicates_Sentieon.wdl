###############################################################################################

#       This WDL script is performs PICARD Mark Duplicates on input sorted BAM Files     ##

###############################################################################################

import "DesignBlock_1/Tasks/PicardMarkDuplicates_Sentieon.wdl" as PICARDMARKDUPLICATES

workflow CallMarkDuplicatesTask {

   call MarkDuplicatesTask

}

