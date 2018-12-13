###########################################################################################

##              This WDL script merges VCFs produced by strelka and mutect               ##

############################################################################################

task combineVariantsTask  {

   File StrelkaVcfBgz                             # Input strelka VCF
   File StrelkaVcfBgzTbi                          # Input strelka VCF index
   File MutectVcfBgz                              # Input mutect VCF
   File MutectVcfBgzTbi                           # Input mutect VCF index

   File Ref                                       # Reference genome fasta
   File RefFai                                    # Reference Genome index

   String SampleName                              # Name of the Sample

   String GatkJar                                 # Path to GATK .jar
   String Java                                    # Path to Java to invoke GATK
   String CombineVariantsThreads                  # No of Threads for the Tool
   String CombineVariantsExtraOptionsString       # String of extra options for the tool, can be empty string

   String Bgzip                                   # Path to BGZip executable

   File BashPreamble                              # bash script to source before every task
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File CombineVariantsScript                     # Path to bash script called within WDL script
   File CombineVariantsEnvProfile                 # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode

   String CombineVariantsSoftMemLimit             # Soft memory limit - nice shutdown
   String CombineVariantsHardMemLimit             # Hard memory limit - kill immediately



   command <<<
        source ${BashPreamble}
        /bin/bash ${CombineVariantsScript} -s ${SampleName} -S ${StrelkaVcfBgz} -M ${MutectVcfBgz} -g ${Ref} -G ${GatkJar} -J ${Java} -Z ${Bgzip} -t ${CombineVariantsThreads} -F ${BashSharedFunctions} -e ${CombineVariantsEnvProfile} -o "'${CombineVariantsExtraOptionsString}'" ${DebugMode}
   >>>


   runtime {
      cpu: "${CombineVariantsThreads}"
      s_vmem: "${CombineVariantsSoftMemLimit}"
      h_vmem: "${CombineVariantsHardMemLimit}"
   }


  output {
      File OutputVcf = "somaticvariants.vcf.gz"
      File OutputVcfIdx = "somaticvariants.vcf.gz.idx"
   }

}
