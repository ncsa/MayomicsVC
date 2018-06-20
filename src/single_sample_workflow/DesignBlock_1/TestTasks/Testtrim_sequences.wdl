##############################################################

###       This WDL script calls the CutAdapt WDL Task       ##

##############################################################

import "DesignBlock_1/Tasks/trim_sequences.wdl" as TRIMSEQ

workflow RunTrimInputSequencesTask {

   call TRIMSEQ.trimsequencesTask

} 
