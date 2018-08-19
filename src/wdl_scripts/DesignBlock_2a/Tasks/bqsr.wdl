###########################################################################################

##                     This WDL script performs BQSR using Sentieon                     ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -b        "Input Deduped Bam"                                     (Required)
#       -k        "List of Known Sites"                                   (Required)
#       -D        "Path to the DBSNP File"                                (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -L        "Sentieon License File"                                 (Required)

############################################################################################

task bqsrTask {

   File InputAlignedSortedDedupedRealignedBam           # Input Sorted Deduped Bam
   File InputAlignedSortedDedupedRealignedBamIdx        # Input Sorted Deduped Bam Index
   File Ref                                             # Reference Genome

   File RefFai                                          # Reference files that are provided as implicit inputs
                                                        # to the WDL Tool to help perform the realignment

   String SampleName                                    # Name of the Sample

   File KnownSites                                      # List of known sites
   File KnownSitesIdx                                   # Index file for the known sites

   File DBSNP                                           # DBSNP file
   File DBSNPIdx                                        # Index file for the DBSNPs   

   String SentieonLicense                               # Sentieon License Information
   String Sentieon                                      # Path to Sentieon
   String SentieonThreads                               # No of Threads for the Tool

   String DebugMode                                     # Enable or Disable Debug Mode

   File BqsrScript                                      # Path to bash script called within WDL script

   command {

      /bin/bash ${BqsrScript} -s ${SampleName}  -L ${SentieonLicense} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -b ${InputAlignedSortedDedupedRealignedBam} -D ${DBSNP} -k ${KnownSites} ${DebugMode}

   }

   
   output {

      File RecalTable = "${SampleName}.recal_data.table"

   }

}
