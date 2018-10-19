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
#       -e        "Path to the environmental profile               (Required)
#       -d        "debug mode on/off                               (Optional: can be empty)


############################################################################################

task bqsrTask {

   File InputAlignedSortedDedupedRealignedBam           # Input Sorted Deduped Bam
   File InputAlignedSortedDedupedRealignedBamBai        # Input Sorted Deduped Bam Index
   File Ref                                             # Reference Genome

   File RefFai                                          # Reference files that are provided as implicit inputs
                                                        # to the WDL Tool to help perform the realignment

   String SampleName                                    # Name of the Sample

   String BQSRKnownSites                                # List of known sites, including dbSNP

   String Sentieon                                      # Path to Sentieon
   String SentieonThreads                               # No of Threads for the Tool

   File BqsrScript                                      # Path to bash script called within WDL script
   File BqsrEnvProfile                                  # File containing the environmental profile variables

   String DebugMode                                     # Enable or Disable Debug Mode


   command {
      /bin/bash ${BqsrScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -b ${InputAlignedSortedDedupedRealignedBam} -k ${BQSRKnownSites} -e ${BqsrEnvProfile} ${DebugMode}
   }

   
   output {
      File RecalTable = "${SampleName}.recal_data.table"
   }

}
