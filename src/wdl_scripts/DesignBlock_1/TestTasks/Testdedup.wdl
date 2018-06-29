#########################################################################

#       This WDL script calls the Marking Duplicates WDL Task     ##

#########################################################################

import "src/wdl_scripts/DesignBlock_1/Tasks/dedup.wdl" as DEDUP

workflow CalldedupTask {

   call DEDUP.dedupTask

}

