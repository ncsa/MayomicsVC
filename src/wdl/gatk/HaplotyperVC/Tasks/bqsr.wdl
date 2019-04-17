###########################################################################################

##                     This WDL script performs BQSR using Sentieon                      ##

##                                Script Options
#       -t        "Number of Threads"                              (Optional)
#       -G        "Reference Genome"                               (Required)
#       -b        "Input Deduped Bam"                              (Required)
#       -k        "List of Known Sites"                            (Required)
#       -D        "Path to the DBSNP File"                         (Required)
#       -s        "Name of the sample"                             (Optional)
#       -S        "Path to the Sentieon Tool"                      (Required)
#       -e        "Path to the environmental profile"              (Required)
#	-F        "Path to the shared functions file"              (Required)
#       -d        "debug mode on/off"                              (Optional: can be empty)


############################################################################################

task bqsrTask {
   String SampleName                     # Name of the Sample

   File InputBams                        # Input Sorted Deduped Bam
   File InputBais                        # Input Sorted Deduped Bam Index

   File Ref                              # Reference Genome
   File RefFai                           # Reference files that are provided as implicit inputs
   File RefDict                          # to the WDL Tool to help perform BQSR


   String BqsrKnownSites                 # List of known sites, including dbSNP
   String GenomicInterval                # Array of chromosome names or genomic intervals for parallel analysis
   File GATKExe                          # Path to GATK4 executable
   String ApplyBQSRExtraOptionsString    # String of extra options for ApplyBQSR. This can be an empty string
   File JavaExe                          # Path to Java8 executable
   String JavaOptionsString              # String of java vm options, like garbage collection and maximum and minimum memory. Can NOT be empty

   String BqsrSoftMemLimit               # Soft memory limit - nice shutdown
   String BqsrHardMemLimit               # Hard memory limit - kill immediately

   File BashPreamble                     # Bash script to run before every task
   File BashSharedFunctions              # Bash script with shared functions
   File BqsrScript                       # Path to bash script called within WDL script

   String DebugMode                      # Enable or Disable Debug Mode

   command <<<
       source ${BashPreamble}
       /bin/bash ${BqsrScript} -s ${SampleName} -b ${InputBams} -G ${Ref} -k ${BqsrKnownSites} -I ${GenomicInterval} -S ${GATKExe} -o ${ApplyBQSRExtraOptionsString} -J ${JavaExe} -e "'${JavaOptionsString}'" -F ${BashSharedFunctions} ${DebugMode}
   >>>
   
   runtime {
      s_vmem: "${BqsrSoftMemLimit}"
      memory: "${BqsrHardMemLimit}"
      docker : "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
   }

   output {
      File OutputBams = "${SampleName}.${GenomicInterval}.bam"
      File OutputBais = "${SampleName}.${GenomicInterval}.bai"
   }
}
