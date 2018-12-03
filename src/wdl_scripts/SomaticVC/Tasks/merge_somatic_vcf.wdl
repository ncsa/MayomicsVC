###########################################################################################

##              This WDL script merges VCFs produced by strelka and mutect               ##

############################################################################################

task mergeSomaticVcfTask  {

   File InputStrelkaVcf                           # Input strelka VCF
   File InputStrelkaVcfIdx                        # Input strelka VCF index
   File InputMutectVcf                            # Input mutect VCF
   File InputMutectVcfIdx                         # Input mutect VCF index

   String SampleName                              # Name of the Sample

   File BashPreamble                              # bash script to source before every task
   File MergeSomaticVcfScript                     # Path to bash script called within WDL script
   File MergeSomaticVcfEnvProfile                 # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${MergeSomaticVcfScript} -s ${SampleName} -S ${InputStrelkaVcf} -M ${InputMutectVcf} -e ${MergeSomaticVcfEnvProfile} ${DebugMode}
   >>>

  output {
      File OutputVcf = "${SampleName}.vcf"
      File OutputVcfIdx = "${SampleName}.vcf.idx"
   }

}
