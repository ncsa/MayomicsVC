###########################################################################################

##                     This WDL script performs VQSR using Sentieon                      ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -V        "Input VCF"                                             (Required)
#       -H        "HapMap Known Sites"                                    (Required)
#       -O        "Omni Known Sites"                                      (Required)
#       -T        "1000 genomes Known Sites"                              (Required)
#       -m        "Mills Known Sites"                                     (Required)
#       -D        "Path to the DBSNP File"                                (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -L        "Sentieon License File"                                 (Required)
#       -d        "debug mode"                                            (Optional)





############################################################################################

task vqsrTask {

   File InputVCF                                         # Input VCF from Haplotyper
   File InputVCFIdx                                      # Input VCF Index

   File Ref                                              # Reference Genome
   File RefFai                                           # Reference files that are provided as implicit inputs

   String SampleName                                     # Name of the Sample

   File HapMapVCF                                        # HapMap known sites
   File HapMapVCFIdx                                     # Index file for HapMap

   File OmniVCF                                          # Omni known sites
   File OmniVCFIdx                                       # Index file for Omni

   File ThousandGVCF                                     # Omni known sites
   File ThousandGVCFIdx                                  # Index file for Omni

   File MillsVCF                                         # Omni known sites
   File MillsVCFIdx                                      # Index file for Omni

   File DBSNP                                            # DBSNP file
   File DBSNPIdx                                         # Index file for the DBSNP   

   String SentieonThreads                                # No of Threads for the Tool
   String SentieonLicense                                # Sentieon License Information
   String Sentieon                                       # Path to Sentieon
   String DebugMode                                      # Enable or Disable Debug Mode

   File VqsrScript                                       # Path to bash script called within WDL script


   command {

      /bin/bash ${VqsrScript} -s ${SampleName}  -L ${SentieonLicense} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -V ${InputVCF} -H ${HapMapVCF} -O ${OmniVCF} -T ${ThousandGVCF} -D ${DBSNP} -m ${MillsVCF} ${DebugMode}

   }

   
   output {
      File RecalibratedVCF = "${SampleName}.recalibrated.vcf"
   }

}
