#################################################################################################

##              This WDL script delivers the results of SomaticVC block                        ##

#################################################################################################

task deliverSomaticVCTask {

   File InputVcf                          # VCF File
   File InputVcfIdx                       # VCF.IDX File

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File BashPreamble                      # Bash script that helps control zombie processes
   File BashSharedFunctions               # Bash script that contains shared helpful functions
   File DeliverySomaticVC_Script          # Bash script that performs the delivery

   String DeliveryFolder_SomaticVC        # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on; optional

   command {
      source ${BashPreamble}
      /bin/bash ${DeliverySomaticVC_Script} -s ${SampleName} -r ${InputVcf} -j ${WorkflowJson} -f ${DeliveryFolder_SomaticVC} -F ${BashSharedFunctions} ${DebugMode}
   }

}
