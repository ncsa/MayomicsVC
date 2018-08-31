###########################################################################################

##                     This WDL script performs VQSR using Sentieon                      ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -V        "Input VCF"                                             (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -e        "Path to the environmental profile                      (Required)
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
   String AnnotateText                  # Annotation text string for vqr 

   String SentieonThreads               # No of Threads for the Tool
   String Sentieon                      # Path to Sentieon

   File VqsrScript                      # Path to bash script called within WDL script
   File VqsrEnvProfile                  # File containing the environmental profile variables

   String DebugMode                     # Enable or Disable Debug Mode


   command {
      /bin/bash ${VqsrScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -V ${InputVCF} -r ${VqsrSnpResourceString} -R ${VqsrIndelResourceString} -a ${AnnotateText} -e ${VqsrEnvProfile} ${DebugMode}
   }

   
   output {
      File RecalibratedVcf = "${SampleName}.INDEL.SNP.recaled.vcf"
      File RecalibratedVcfIdx = "${SampleName}.INDEL.SNP.recaled.vcf.idx"
   }

}
