#################################################################################################

##              This WDL script delivers output of alignment Block                             ##

#                              Script Options
#      -b        "BAM File"                                  (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverAlignmentTask {

   File InputBams                         # aligned sorted dedupped BAM file
   File InputBais                         # aligned sorted dedupped BAM.BAI file

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File DeliveryAlignment_Script          # Bash script that performs the delivery
   String DeliveryFolder_Alignment        # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      /bin/bash ${DeliveryAlignment_Script} -s ${SampleName} -b ${AlignedSortedDedupedBam} -j ${WorkflowJson} -f ${DeliveryFolder_Alignment} ${DebugMode}
   }

}
