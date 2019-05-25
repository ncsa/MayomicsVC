#################################################################################################
#                                                                                               #  
##              This WDL script merges the input asligned BAMs before deduplication            ##
#                                                                                               #  
#################################################################################################

task mergebamTask {
   String SampleName               # Name of the Sample

   Array[File] InputBams           # Input Sorted BAM File
   Array[File] InputBais           # Input Sorted Bam Index File


   String SamtoolsExe              # Path to Samtools Executable

   String DebugMode                # Variable to check whether Debud Mode is on

   String MergeSoftMemLimit        # Soft memory limit - nice shutdown
   String MergeHardMemLimit        # Hard memory limit - kill immediately
   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File MergeBamScript             # Bash script that is called inside the WDL script

   command <<<
   	   #source ${BashPreamble}
	   module load BWA/0.7.17-IGB-gcc-4.9.4 SAMtools/1.7-IGB-gcc-4.9.4 GATK/4.0.9.0-IGB-gcc-4.9.4-Java-1.8.0_152-Python-3.6.1
   	   /bin/bash ${MergeBamScript} -b ${sep=',' InputBams} -s ${SampleName} -S ${SamtoolsExe} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      s_vmem: "${MergeSoftMemLimit}"
      memory: "${MergeHardMemLimit}"
      docker : "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }
}
