##############################################################

###       This WDL script calls the CutAdapt WDL Task       ##

##############################################################

import "src/wdl_scripts/Alignment/Tasks/trim_sequences.wdl" as TRIMSEQ

workflow RunTrimInputSequencesTask {

   call TRIMSEQ.trimsequencesTask

} 
