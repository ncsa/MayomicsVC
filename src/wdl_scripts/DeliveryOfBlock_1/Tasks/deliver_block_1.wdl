#################################################################################################

##              This WDL script delivers output of Design Block 1                              ##

#                              Script Options
#      -b        "BAM File"                                  (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverBlock1Task {

   File AlignedSortedDedupedBam           # aligned sorted dedupped BAM file
   File AlignedSortedDedupedBamBai        # aligned sorted dedupped BAM.BAI file

   String SampleName                      # Name of the Sample

   File DeliveryBlock_1_Script            # Bash script that performs the delivery
   String DeliveryFolder_Block_1          # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      /bin/bash ${DeliveryBlock_1_Script} -s ${SampleName} -b ${AlignedSortedDedupedBam} -f ${DeliveryFolder_Block_1} ${DebugMode}
   }

}
