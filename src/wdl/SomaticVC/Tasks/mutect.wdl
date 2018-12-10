##########################################################################################

##              This WDL script performs tumor/normal Variant Calling  using mutect     ##

##########################################################################################

task mutectTask {

   File TumorBams                                 # Input Sorted Deduped Tumor Bam
   File TumorBais                                 # Input Sorted Deduped Tumor Bam Index
   File NormalBams                                # Input Sorted Deduped Normal Bam
   File NormalBais                                # Input Sorted Deduped Normal Bam Index

   File Ref                                       # Reference Genome
   File RefFai                                    # Reference Genome index
   File RefDict                                   # Reference Dictionary

   String SampleName                              # Name of the Sample

   String MutectExtraOptionsString                # String of extra options for mutect, this can be an empty string

   String GatkJarPath                             # Path to GATK .jar
   String JavaPath                                # Path to Java to invoke GATK
   String MutectThreads                           # No of Threads for the Tool
   String MutectJavaMemOption                     # is a string in form of e.g. "-Xmx4G" or "-Xms2G -Xmx8G"

   String BcfToolsPath                            # Path to BCF tools executable
   String BgzipPath                               # Path to BGZip executable
   String SamtoolsPath                            # Path to Samtools executable

   File BashPreamble                              # Bash script that helps control zombie processes
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File MutectScript                              # Path to bash script called within WDL script
   File MutectEnvProfile                          # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode

   String MutectSoftMemLimit                     # Soft memory limit - nice shutdown
   String MutectHardMemLimit                     # Hard memory limit - kill immediately


   command <<<
        source ${BashPreamble}
        /bin/bash ${MutectScript} -s ${SampleName} -G ${GatkJarPath} -J ${JavaPath} -j ${MutectJavaMemOption} -B ${BcfToolsPath} -Z ${BgzipPath} -S ${SamtoolsPath} -g ${Ref} -t ${MutectThreads} -T ${TumorBams} -N ${NormalBams} -o ${MutectExtraOptionsString} -e ${MutectEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>


   runtime {
      cpu: "${MutectThreads}"
      s_vmem: "${MutectSoftMemLimit}"
      h_vmem: "${MutectHardMemLimit}"
   }


  output {
      File OutputVcf = "${SampleName}.vcf"
      File OutputVcfIdx = "${SampleName}.vcf.idx"
   }

}
