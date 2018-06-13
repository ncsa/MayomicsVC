###############################################################################################

#       This WDL script is performs PICARD Mark Duplicates on input sorted BAM Files     ##

###############################################################################################

import "DesignBlock_1/Tasks/Sentieon/dedup.wdl" as DEDUP

workflow CalldedupTask {

   call DEDUP.dedupTask

}

