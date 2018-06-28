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
#               -L        "Sentieon License File"                                 (Required)

############################################################################################

task variantCallingTask {

   File InputAlignedSortedDedupedRealignedBam              # Input Sorted Deduped Bam
   File InputAlignedSortedDedupedRealignedBamIdx           # Input Sorted Deduped Bam Index

   File Ref                                                # Reference Genome

   String SampleName                                       # Name of the Sample

   File DBSNP                                              # DBSNP file
   File DBSNPIdx                                           # Index file for the DBSNPs   
  
   File RecalTable                                         # Recal Table after BQSR step
   
   String Threads                                          # No of Threads for the Tool
   String SentieonLicense                                  # Sentieon License Information
   String Sentieon                                         # Path to Sentieon
   String DebugMode                                        # Enable or Disable Debug Mode

   File HaplotyperScript                                   # Path to bash script called within WDL script

   command {

      /bin/bash ${HaplotyperScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${Threads} -b ${InputAlignedSortedDedupedRealignedBam} -D ${DBSNP} -r ${RecalTable} -L ${SentieonLicense} ${DebugMode}

   }

  output {
   
      File VCF = "${SampleName}.vcf"
    
   }

}
