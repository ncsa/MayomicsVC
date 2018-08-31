#################################################################################################

##              This WDL script delivers the results of Design block 2a                        ##

#                              Script Options
#      -s        "Input SNP VCF File"                        (Required)      
#      -i        "Input Indel VCF File"                      (Required)      
#      -j        "Input JSON File"                           (Required)      
#      -f        "Path to the delivery folder "              (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task deliverBlock2aTask {

   File RecalibratedSnpVcf                # Input SNP VCF File
   File RecalibratedSnpVcfIdx             # Input SNP VCF.IDX File
   File RecalibratedIndelVcf              # Input Indel VCF File
   File RecalibratedIndelVcfIdx           # Input Indel VCF.IDX File

   File InputJson                         # Input JSON file for the workflow

   File DeliveryBlock_2a_Script           # Bash script that performs the delivery
   File DeliveryFolder_Block_2a           # Input delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      /bin/bash ${DeliveryBlock_2a_Script} -s ${RecalibratedSnpVcf} -i ${RecalibratedIndelVcf} -j ${InputJson} -f ${DeliveryFolder_Block_2a} ${DebugMode}
   }

}
