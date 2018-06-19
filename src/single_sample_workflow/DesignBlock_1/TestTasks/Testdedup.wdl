#########################################################################

#       This WDL script calls the Marking Duplicates WDL Task     ##

#########################################################################

import "DesignBlock_1/Tasks/dedup.wdl" as DEDUP

workflow CalldedupTask {

   call DEDUP.dedupTask

}

