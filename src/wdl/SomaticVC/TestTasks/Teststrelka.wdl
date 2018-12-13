#########################################################################                                                          
#       This WDL script calls the Strelka WDL Task                     ##

#########################################################################

import "MayomicsVC/src/wdl/SomaticVC/Tasks/strelka.wdl" as STRELKA

workflow CallStrelkaTask {

   call STRELKA.strelkaTask as strelka

   output {
      File StrelkaVcf = strelka.OutputVcfBgz
      File StrelkaVcfIdx = strelka.OutputVcfBgzTbi
   }

}

