#################################################################################################
#                                                                                               #  
##       This WDL script merges the input gvcfs  before joint calling multiple samples         ##
#                                                                                               #  
#################################################################################################

task mergegvcfsTask {
   String SampleName                              # Name of the Sample

   Array[File] InputGvcfs                         # Input GVCF files
   Array[File] InputIdxs                          # Input GVCF files Index

   File GATKExe                                   # Path to GATK4 executable
   File JavaExe                                   # Path to Java8 executable
   String JavaOptionsString                       # String of java vm options. Can NOT be empty

   String MergeSoftMemLimit                       # Soft memory limit - nice shutdown
   String MergeHardMemLimit                       # Hard memory limit - kill immediately
   File BashPreamble                              # Bash script that helps control zombie processes
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File MergeGvcfsScript                          # Bash script that is called inside the WDL script

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        #source ${BashPreamble}
	module load BWA/0.7.17-IGB-gcc-4.9.4 SAMtools/1.7-IGB-gcc-4.9.4 GATK/4.0.9.0-IGB-gcc-4.9.4-Java-1.8.0_152-Python-3.6.1
        /bin/bash ${MergeGvcfsScript} -s ${SampleName} -b ${sep=',' InputGvcfs} -S ${GATKExe} -J ${JavaExe} -e "'${JavaOptionsString}'" -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      s_vmem: "${MergeSoftMemLimit}"
      memory: "${MergeHardMemLimit}"
      docker : "broadinstitute/gatk:4.1.1.0"
   }

  output {
      File OutputVcf = "${SampleName}.g.vcf"
      File OutputVcfIdx = "${SampleName}.g.vcf.idx"
   }

}
