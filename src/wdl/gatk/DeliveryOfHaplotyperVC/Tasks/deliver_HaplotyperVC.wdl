#################################################################################################

##              This WDL script delivers the results of HaplotyperVC block                     ##

#################################################################################################

task deliverHaplotyperVCTask {

   File InputVcf                          # VCF File
   File InputVcfIdx                       # VCF.IDX File

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File BashPreamble                      # Bash script that helps control zombie processes
   File BashSharedFunctions               # Bash script that contains shared helpful functions
   File DeliveryHaplotyperVC_Script       # Bash script that performs the delivery

   String DeliveryFolder_HaplotyperVC     # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on; optional

   command {
      source ${BashPreamble}
      /bin/bash ${DeliveryHaplotyperVC_Script} -s ${SampleName} -r ${InputVcf} -j ${WorkflowJson} -f ${DeliveryFolder_HaplotyperVC} -F ${BashSharedFunctions} ${DebugMode}
   }

}
