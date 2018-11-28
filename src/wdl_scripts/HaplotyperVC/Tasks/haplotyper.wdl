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

   File InputBams                                 # Input Sorted Deduped Bam
   File InputBais                                 # Input Sorted Deduped Bam Index
   File RecalTable                                # Input Recal Table after BQSR step

   File Ref                                       # Reference Genome
   File RefFai                                    # Reference Genome index

   String SampleName                              # Name of the Sample

   String HaplotyperExtraOptionsString            # String of extra options for haplotyper, this can be an empty string

   File DBSNP                                     # DBSNP file
   File DBSNPIdx                                  # Index file for the DBSNPs


   String Sentieon                                # Path to Sentieon
   String SentieonThreads                         # No of Threads for the Tool

   String HaplotyperSoftMemLimit                  # Soft memory limit - nice shutdown
   String HaplotyperHardMemLimit                  # Hard memory limit - kill immediately

   File BashPreamble                              # bash script to source before every task
   File BashSharedFunctions                       # Bash script with shared functions
   File HaplotyperScript                          # Path to bash script called within WDL script
   File HaplotyperEnvProfile                      # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${HaplotyperScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -b ${InputBams} -D ${DBSNP} -r ${RecalTable} -o "'${HaplotyperExtraOptionsString}'" -e ${HaplotyperEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${HaplotyperSoftMemLimit}"
      h_vmem: "${HaplotyperHardMemLimit}"
   }

  output {
   
      File OutputVcf = "${SampleName}.vcf"
      File OutputVcfIdx = "${SampleName}.vcf.idx"
    
   }

}
