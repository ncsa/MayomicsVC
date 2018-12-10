#################################################################################################

##              This WDL script delivers output of Alignment Block                             ##

#################################################################################################

task deliverAlignmentTask {

   File InputBams                         # aligned sorted dedupped BAM file
   File InputBais                         # aligned sorted dedupped BAM.BAI file

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File BashPreamble                      # Bash script that helps control zombie processes
   File BashSharedFunctions               # Bash script that contains shared helpful functions
   File DeliveryAlignment_Script          # Bash script that performs the delivery

   String DeliveryFolder_Alignment        # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      source ${BashPreamble}
      /bin/bash ${DeliveryAlignment_Script} -s ${SampleName} -b ${InputBams} -j ${WorkflowJson} -f ${DeliveryFolder_Alignment} -F ${BashSharedFunctions} ${DebugMode}
   }

}
