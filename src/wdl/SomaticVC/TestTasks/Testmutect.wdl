#########################################################################                                                          
#       This WDL script calls the Mutect WDL Task                      ##

#########################################################################

import "src/wdl/SomaticVC/Tasks/mutect.wdl" as STRELKA

workflow CallMutectTask {

   call STRELKA.mutectTask as mutect

   output {
      Array[File] MutectVcf = mutect.OutputVcfBgz
      Array[File] MutectVcfIdx = mutect.OutputVcfBgzTbi
   }

}

