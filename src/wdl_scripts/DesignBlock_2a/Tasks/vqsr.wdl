###########################################################################################

##                     This WDL script performs VQSR using Sentieon                      ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -V        "Input VCF"                                             (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -L        "Sentieon License File"                                 (Required)
#       -d        "debug mode"                                            (Optional)





############################################################################################

task vqsrTask {

   File InputVCF                        # Input VCF from Haplotyper
   File InputVCFIdx                     # Input VCF Index
  
   File Ref                             # Reference Genome
   File RefFai                          # Reference files that are provided as implicit inputs

   String SampleName                    # Name of the Sample

   String VqsrSnpResourceString         # Sentieon resource string that contains paths to 
                                        # 1000G, omni, dbSNP, hapmap VCFs and associated parameter values
   String VqsrIndelResourceString       # Sentieon resource string that contains paths to 
                                        # dbSNP, Mills VCFs and associated parameter values

   String SentieonThreads               # No of Threads for the Tool
   String SentieonLicense               # Sentieon License Information
   String Sentieon                      # Path to Sentieon
   String DebugMode                     # Enable or Disable Debug Mode

   File VqsrScript                      # Path to bash script called within WDL script


   command {
#      /bin/bash ${VqsrScript} -s ${SampleName}  -L ${SentieonLicense} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -V ${InputVCF} -r ${VqsrSnpResourceString} -R ${VqsrIndelResourceString} ${DebugMode}
      /bin/bash ${VqsrScript} -s ${SampleName}  -L ${SentieonLicense} -S ${Sentieon} -G ${Ref} -V ${InputVCF} -r ${VqsrSnpResourceString} -R ${VqsrIndelResourceString} ${DebugMode}
   }

   
   output {
      File RecalibratedVCF = "${SampleName}.recalibrated.vcf"
   }

}
