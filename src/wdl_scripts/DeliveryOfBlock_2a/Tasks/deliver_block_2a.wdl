#################################################################################################

##              This WDL script delivers the results of Design block 2a                        ##

#                              Script Options
#      -s        "SNP VCF File"                        (Required)      
#      -i        "Indel VCF File"                      (Required)      
#      -j        "JSON File"                           (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverBlock2aTask {

   File RecalibratedSnpVcf                # SNP VCF File
   File RecalibratedSnpVcfIdx             # SNP VCF.IDX File
   File RecalibratedIndelVcf              # Indel VCF File
   File RecalibratedIndelVcfIdx           # Indel VCF.IDX File

   File InputJson                         # JSON file for the workflow

   File DeliveryBlock_2a_Script           # Bash script that performs the delivery
   String DeliveryFolder_Block_2a         # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      /bin/bash ${DeliveryBlock_2a_Script} -s ${RecalibratedSnpVcf} -i ${RecalibratedIndelVcf} -j ${InputJson} -f ${DeliveryFolder_Block_2a} ${DebugMode}
   }

}
