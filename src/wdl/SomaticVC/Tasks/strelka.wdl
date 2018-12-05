###########################################################################################

##              This WDL script performs tumor/normal Variant Calling  using strelka     ##

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

   File BashPreamble                              # Bash script that helps control zombie processes
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File StrelkaScript                             # Path to bash script called within WDL script
   File StrelkaEnvProfile                         # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${StrelkaScript} -s ${SampleName} -S ${Strelka} -G ${Ref} -t ${StrelkaThreads} -T ${TumorBams} -N ${NormalBams} -o ${StrelkaExtraOptionsString} -e ${StrelkaEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

  output {
      File OutputVcf = "${SampleName}.vcf"
      File OutputVcfIdx = "${SampleName}.vcf.idx"
   }

}
