#################################################################################################

##              This WDL script delivers the results of HaplotyperVC block                        ##

#                              Script Options
#      -s        "Sample name"                               (Required)      
#      -r        "Recalibrated VCF File"                     (Required)      
#      -j        "JSON File"                                 (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverHaplotyperVCTask {

   File InputVcf                          # VCF File
   File InputVcfIdx                       # VCF.IDX File

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File BashPreamble
   File DeliveryHaplotyperVC_Script       # Bash script that performs the delivery
   String DeliveryFolder_HaplotyperVC     # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      source ${BashPreamble}
      /bin/bash ${DeliveryHaplotyperVC_Script} -s ${SampleName} -r ${InputVcf} -j ${WorkflowJson} -f ${DeliveryFolder_HaplotyperVC} ${DebugMode}
   }

}
