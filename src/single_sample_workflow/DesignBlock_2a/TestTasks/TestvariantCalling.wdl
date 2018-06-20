#######################################################################

###       This WDL script calls the Variant calling WDL Task       ##

#######################################################################

import "DesignBlock_2a/Tasks/variantCalling.wdl" as HAPLOTYPER

workflow CallvariantCallingTask {

   call HAPLOTYPER.variantCallingTask

}

