#################################################################################################

##              This WDL script marks the duplicates on input sorted BAMs                ##

#                              Script Options
#      -b        "Input BAM File"                            (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverBlock1Task {

   File AlignedSortedDedupedBam           # Input aligned sorted dedupped BAM file
   File AlignedSortedDedupedBamBai        # Input aligned sorted dedupped BAM.BAI file

   File DeliveryBlock_1_Script            # Bash script that performs the delivery
   File DeliveryFolder_Block_1            # Input delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      /bin/bash ${DeliveryBlock_1_Script} -b ${AlignedSortedDedupedBam} -f ${DeliveryFolder_Block_1} ${DebugMode}
   }

}
