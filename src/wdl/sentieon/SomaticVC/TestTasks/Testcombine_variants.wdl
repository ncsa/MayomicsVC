#########################################################################                                                                   
#       This WDL script calls the Mutect combine_variants Task         ##

#########################################################################

import "MayomicsVC/src/wdl/sentieon/SomaticVC/Tasks/combine_variants.wdl" as COMBINE_VARIANTS

workflow CallCombineVariantsTask {

   call COMBINE_VARIANTS.combineVariantsTask as CombineVariants

   output {
      File OutputVcfGz = CombineVariants.OutputVcfGz
      File OutputVcfGzTbi = CombineVariants.OutputVcfGzTbi
   }

}

