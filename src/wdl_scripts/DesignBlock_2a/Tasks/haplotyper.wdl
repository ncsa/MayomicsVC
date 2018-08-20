###########################################################################################

##              This WDL script performs Variant Calling  using Sentieon                 ##

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
   File InputAlignedSortedDedupedRealignedBamBai           # Input Sorted Deduped Bam Index
   File RecalTable                                         # Input Recal Table after BQSR step

   File Ref                                                # Reference Genome

   String SampleName                                       # Name of the Sample

   File DBSNP                                              # DBSNP file
   File DBSNPIdx                                           # Index file for the DBSNPs   
  
   
   String SentieonLicense                                  # Sentieon License Information
   String Sentieon                                         # Path to Sentieon
   String SentieonThreads                                  # No of Threads for the Tool

   String DebugMode                                        # Enable or Disable Debug Mode

   File HaplotyperScript                                   # Path to bash script called within WDL script

   command {

      /bin/bash ${HaplotyperScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -b ${InputAlignedSortedDedupedRealignedBam} -D ${DBSNP} -r ${RecalTable} -L ${SentieonLicense} ${DebugMode}

   }

  output {
   
      File VCF = "${SampleName}.vcf"
      File VcfIdx = "${SampleName}.vcf.idx"
    
   }

}
