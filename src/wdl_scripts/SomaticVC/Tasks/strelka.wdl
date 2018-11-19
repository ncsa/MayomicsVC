###########################################################################################

##              This WDL script performs tumor/normal Variant Calling  using strelka     ##

##                                Script Options
#               -t        "Number of Threads"                                     (Optional)
#               -G        "Reference Genome"                                      (Required)
#               -T        "Input Sorted Deduped Tumor Bam"                        (Required)
#               -N        "Input Sorted Deduped Normal Bam"                       (Required)
#               -s        "Name of the sample"                                    (Optional)
#               -o        "Strelka Extra Options"                                 (Required)
#               -S        "Path to the Strelka Tool"                              (Required)
#               -e        "Path to the environmental profile                      (Required)
#               -d        "debug mode on/off                        (Optional: can be empty)
#

############################################################################################

task strelkaTask {

   File TumorBams                                 # Input Sorted Deduped Tumor Bam
   File TumorBais                                 # Input Sorted Deduped Tumor Bam Index
   File NormalBams                                # Input Sorted Deduped Normal Bam
   File NormalBais                                # Input Sorted Deduped Normal Bam Index

   File Ref                                       # Reference Genome
   File RefFai                                    # Reference Genome index

   String SampleName                              # Name of the Sample

   String StrelkaExtraOptionsString               # String of extra options for strelka, this can be an empty string

   String Strelka                                 # Path to Strelka 
   String StrelkaThreads                          # No of Threads for the Tool

   File BashPreamble                              # bash script to source before every task
   File StrelkaScript                             # Path to bash script called within WDL script
   File StrelkaEnvProfile                         # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${StrelkaScript} -s ${SampleName} -S ${Sentieon} -G ${Ref} -t ${SentieonThreads} -T ${TumorBams} -N ${NormalBams} -o ${StrelkaExtraOptionsString} -e ${StrelkaEnvProfile} ${DebugMode}
   >>>

  output {
      File OutputVcf = "${SampleName}.vcf"
      File OutputVcfIdx = "${SampleName}.vcf.idx"
   }

}
