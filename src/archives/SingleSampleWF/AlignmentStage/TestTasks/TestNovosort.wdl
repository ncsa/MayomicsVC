#######################################################################################

##       This WDL script is performs Novosort to sort the input BAM Files     ##

#######################################################################################

import "AlignmentStage_WDL/Tasks/Novosort.wdl" as NOVOSORT

workflow CallNovosortTask {

      call NOVOSORT.NovosortTask
}
   
