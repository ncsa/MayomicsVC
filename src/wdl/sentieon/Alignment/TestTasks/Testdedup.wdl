#########################################################################

#       This WDL script calls the Marking Duplicates WDL Task     ##

#########################################################################

import "MayomicsVC/src/wdl/sentieon/Alignment/Tasks/dedup.wdl" as DEDUP

workflow CalldedupTask {

   call DEDUP.dedupTask

}

