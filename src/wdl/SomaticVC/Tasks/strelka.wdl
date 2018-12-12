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
   File IndelGtFixPerlScript                      # path and name of indel_GT_fix_perl_script
   File SnpGtFixPerlScript                        # path and name of snp_GT_fix_perl_script

   String Bcftools                                # Path to BCFtools
   String Samtools                                # Path to Samtools
   String Bgzip                                   # Path to bgzip

   File BashPreamble                              # Bash script that helps control zombie processes
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File StrelkaScript                             # Path to bash script called within WDL script
   File StrelkaEnvProfile                         # File containing the environmental profile variables

   String DebugMode                               # Enable or Disable Debug Mode

   String StrelkaSoftMemLimit                     # Soft memory limit - nice shutdown
   String StrelkaHardMemLimit                     # Hard memory limit - kill immediately


   command <<<
        source ${BashPreamble}
        /bin/bash ${StrelkaScript} -s ${SampleName} -N ${NormalBams} -T ${TumorBams} -g ${Ref} -B ${Bcftools} -I ${Strelka} -S ${Samtools} -Z ${Bgzip} -t ${StrelkaThreads} -e ${StrelkaEnvProfile} -F ${BashSharedFunctions} -i ${IndelGtFixPerlScript} -p ${SnpGtFixPerlScript} -o "'${StrelkaExtraOptionsString}'" ${DebugMode}
   >>>


   runtime {
      cpu: "${StrelkaThreads}"
      s_vmem: "${StrelkaSoftMemLimit}"
      h_vmem: "${StrelkaHardMemLimit}"
   }


  output {
      File OutputVcfBgz = "${SampleName}.vcf.bgz"
      File OutputVcfBgzTbi = "${SampleName}.vcf.bgz.tbi"
   }

}
