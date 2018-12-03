#######################################################################

###       This WDL script calls the Variant calling WDL Task       ##

#######################################################################

import "src/wdl_scripts/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER

workflow CallvariantCallingTask {

   call HAPLOTYPER.variantCallingTask

}

