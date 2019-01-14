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

   String GatkJar                                 # Path to GATK .jar
   String Java                                    # Path to Java to invoke GATK
   String MutectThreads                           # No of Threads for the Tool
   String MutectJavaMemOption                     # is a string in form of e.g. "-Xmx4G" or "-Xms2G -Xmx8G"
   String MutectExtraOptionsString                # String of extra options for mutect, this can be an empty string

   String Bcftools                                # Path to BCF tools executable
   String Bgzip                                   # Path to BGZip executable
   String Samtools                                # Path to Samtools executable

   File BashPreamble                              # Bash script that helps control zombie processes
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File MutectScript                              # Path to bash script called within WDL script
   File MutectEnvProfile                          # File containing the environmental profile variables
   File FixDPScript                               # Path to fix DP script

   String DebugMode                               # Enable or Disable Debug Mode

   String MutectSoftMemLimit                     # Soft memory limit - nice shutdown
   String MutectHardMemLimit                     # Hard memory limit - kill immediately


   command <<<
        source ${BashPreamble}
        /bin/bash ${MutectScript} -s ${SampleName} -N ${NormalBams} -T ${TumorBams} -g ${Ref} -G ${GatkJar} -J ${Java} -j "'${MutectJavaMemOption}'" -B ${Bcftools} -Z ${Bgzip} -S ${Samtools} -t ${MutectThreads} -F ${BashSharedFunctions} -e ${MutectEnvProfile} -D ${FixDPScript} -o "'${MutectExtraOptionsString}'" ${DebugMode}
   >>>


   runtime {
      cpu: "${MutectThreads}"
      s_vmem: "${MutectSoftMemLimit}"
      h_vmem: "${MutectHardMemLimit}"
   }


  output {
      File OutputVcfBgz = "${SampleName}.vcf.bgz"
      File OutputVcfBgzTbi = "${SampleName}.vcf.bgz.tbi"
  }


}
