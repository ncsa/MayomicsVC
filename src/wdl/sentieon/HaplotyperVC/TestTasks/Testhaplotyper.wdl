#######################################################################

###       This WDL script calls the Variant calling WDL Task       ##

#######################################################################

import "MayomicsVC/src/wdl/sentieon/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER

workflow CallvariantCallingTask {

   call HAPLOTYPER.variantCallingTask

}

