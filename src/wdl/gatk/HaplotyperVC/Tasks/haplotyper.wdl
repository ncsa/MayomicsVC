###########################################################################################

##              This WDL script performs Variant Calling  using Sentieon                 ##

##                                Script Options
#               -t        "Number of Threads"                                     (Optional)
#               -G        "Reference Genome"                                      (Required)
#               -b        "Input Sorted Deduped Bam"                              (Required)
#               -D        "DBSNP File"                                            (Required)
#               -s        "Name of the sample"                                    (Optional)
#               -r        "Recal Data Table"                                      (Required)
#               -o        "Haplotyper Extra Options"                              (Required)
#               -S        "Path to the Sentieon Tool"                             (Required)
#               -e        "Path to the environmental profile                      (Required)
#               -F        "Path to shared functions file"                         (Required)
#               -d        "debug mode on/off                        (Optional: can be empty)
#

############################################################################################

task variantCallingTask {
   String SampleName                              # Name of the Sample

   File InputBams                                 # Input Sorted Deduped Bam
   File InputBais                                 # Input Sorted Deduped Bam Index

   File Ref                                       # Reference Genome
   File RefFai                                    # Reference Genome index
   File RefDict                                   # Reference Genome dictionary

   File DBSNP                                     # DBSNP file
   File DBSNPIdx                                  # Index file for the DBSNPs
   String GenomicInterval                         # Array of chromosome names or genomic intervals for parallel analysis

   File GATKExe                                   # Path to GATK4 executable
   String HaplotyperThreads
   String HaplotyperExtraOptionsString            # String of extra options for haplotyper, this can be an empty string
   File JavaExe                                   # Path to Java8 executable
   String JavaOptionsString                       # String of java vm options. Can NOT be empty

   String HaplotyperSoftMemLimit                  # Soft memory limit - nice shutdown
   String HaplotyperHardMemLimit                  # Hard memory limit - kill immediately

   File BashPreamble                              # bash script to source before every task
   File BashSharedFunctions                       # Bash script with shared functions
   File HaplotyperScript                          # Path to bash script called within WDL script

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${HaplotyperScript} -s ${SampleName} -b ${InputBams} -G ${Ref} -D ${DBSNP} -I ${GenomicInterval} -S ${GATKExe} -t ${HaplotyperThreads} -o "'${HaplotyperExtraOptionsString}'" -J ${JavaExe} -e "'${JavaOptionsString}'" -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${HaplotyperThreads}"
      s_vmem: "${HaplotyperSoftMemLimit}"
      memory: "${HaplotyperHardMemLimit}"
      docker : "broadinstitute/gatk:4.1.1.0"
   }

  output {
      File OutputVcf = "${SampleName}.${GenomicInterval}.g.vcf"
      File OutputVcfIdx = "${SampleName}.${GenomicInterval}.g.vcf.idx"
   }

}
