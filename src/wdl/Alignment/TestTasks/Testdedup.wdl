#########################################################################

#       This WDL script calls the Marking Duplicates WDL Task     ##

#########################################################################

import "src/wdl/Alignment/Tasks/dedup.wdl" as DEDUP

workflow CalldedupTask {

   call DEDUP.dedupTask

}

