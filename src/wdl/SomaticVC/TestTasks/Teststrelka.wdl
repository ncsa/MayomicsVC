#########################################################################                                                          
#       This WDL script calls the Strelka WDL Task                     ##

#########################################################################

import "src/wdl/SomaticVC/Tasks/strelka.wdl" as STRELKA

workflow CallStrelkaTask {

   call STRELKA.strelkaTask as strelka

   output {
      Array[File] StrelkaVcf = strelka.OutputVcfBgz
      Array[File] StrelkaVcfIdx = strelka.OutputVcfBgzTbi
   }

}

