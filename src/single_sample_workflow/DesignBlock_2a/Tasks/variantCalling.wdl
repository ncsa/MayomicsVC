###########################################################################################

##              This WDL script performs Variant Calling  using Sentieon                ##

##                                Script Options
#               -t        "Number of Threads"                                     (Optional)
#               -G        "Reference Genome"                                      (Required)
#               -b        "Input Sorted Deduped Bam"                              (Required)
#               -D        "DBSNP File"                                            (Required)
#               -s        "Name of the sample"                                    (Optional)
#               -r        "Recal Data Table"                                      (Required)
#               -S        "Path to the Sentieon Tool"                             (Required)
#               -O        "Directory for the Output"                              (Required)
#               -L        "Sentieon License File"                                 (Required)
#               -e        "Path to the Error Log File"                            (Required)
#               -d        "Debug Mode Toggle"                                     (Optional)

############################################################################################

task variantCallingTask {

   File InputBam                   # Input Sorted Deduped Bam
   File InputBamIdx                # Input Sorted Deduped Bam Index
   File RefFasta                   # Reference Genome

   String sampleName               # Name of the Sample

   File DBSNP                      # DBSNP file
   File DBSNP_Idx                  # Index file for the DBSNPs   
   
   File RecalTable                 # Recal Table after BQSR step
   
   String Threads                  # No of Threads for the Tool
   String Sentieon_License         # Sentieon License Information
   String Sentieon                 # Path to Sentieon
   String OutDir                   # Output Directory
   Boolean Debug_Mode_EN           # Enable or Disable Debug Mode
   String Error_Logs               # Filepath to Error Logs

   File VariantCalling_Script      # Path to bash script called within WDL script

   command {

      /bin/bash ${VariantCalling_Script} -s ${sampleName} -O ${OutDir} -S ${Sentieon} -G ${RefFasta} -t ${Threads} -b ${InputBam} -D ${DBSNP} -r ${RecalTable} -e ${Error_Logs} -d ${Debug_Mode_EN} -L ${Sentieon_License}

   }


   output {
   
      File VCFFile = "${OutDir}/${sampleName}.vcf"
    
   }

}
