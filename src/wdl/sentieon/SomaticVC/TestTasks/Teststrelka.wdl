#########################################################################                                                          
#       This WDL script calls the Strelka WDL Task                     ##

#########################################################################

import "MayomicsVC/src/wdl/sentieon/SomaticVC/Tasks/strelka.wdl" as STRELKA

workflow CallStrelkaTask {

   call STRELKA.strelkaTask as strelka

   output {
      File StrelkaVcfBgz = strelka.OutputVcfBgz
      File StrelkaVcfBgzTbi = strelka.OutputVcfBgzTbi
   }

}

