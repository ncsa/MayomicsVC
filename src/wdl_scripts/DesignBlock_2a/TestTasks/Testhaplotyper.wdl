#######################################################################

###       This WDL script calls the Variant calling WDL Task       ##

#######################################################################

import "wdl_scripts/DesignBlock_2a/Tasks/haplotyper.wdl" as HAPLOTYPER

workflow CallvariantCallingTask {

   call HAPLOTYPER.variantCallingTask

}

