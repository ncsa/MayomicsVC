#################################################################################################

##              This WDL script delivers the results of Design block 2a                        ##

#                              Script Options
#      -s        "Sample name"                               (Required)      
#      -r        "Recalibrated VCF File"                     (Required)      
#      -j        "JSON File"                                 (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverBlock2aTask {

   File RecalibratedVcf                   # VCF File
   File RecalibratedVcfIdx                # VCF.IDX File

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File DeliveryBlock_2a_Script           # Bash script that performs the delivery
   String DeliveryFolder_Block_2a         # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      /bin/bash ${DeliveryBlock_2a_Script} -s ${SampleName} -r ${RecalibratedVcf} -j ${WorkflowJson} -f ${DeliveryFolder_Block_2a} ${DebugMode}
   }

}
