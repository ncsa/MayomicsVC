###########################################################################################

##                     This WDL script performs VQSR using Sentieon                      ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -V        "Input VCF"                                             (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -e        "Path to the environmental profile                      (Required)
#       -F        "Path to shared functions file"                         (Required)
#       -d        "debug mode"                                            (Optional)
############################################################################################

task vqsrTask {

   File InputVcf                        # Input VCF from Haplotyper
   File InputVcfIdx                     # Input VCF Index
  
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

   String VqsrSoftMemLimit              # Soft memory limit - nice shutdown
   String VqsrHardMemLimit              # Hard memory limit - kill immediately

   File BashPreamble                    # Path to bash script to source before every task
   File BashSharedFunctions             # Bash script with shared functions
   File VqsrScript                      # Path to bash script called within WDL script
   File VqsrEnvProfile                  # File containing the environmental profile variables

   String DebugMode                     # Enable or Disable Debug Mode

   command <<<
      source ${BashPreamble}
      /bin/bash ${VqsrScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -V ${InputVcf} -r "'${VqsrSnpResourceString}'" -R "'${VqsrIndelResourceString}'" -a ${AnnotateText} -e ${VqsrEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${VqsrSoftMemLimit}"
      h_vmem: "${VqsrHardMemLimit}"
   }

   output {
      File OutputVcf = "${SampleName}.SNP.recaled.vcf"
      File OutputVcfIdx = "${SampleName}.SNP.recaled.vcf.idx"
   }
}
