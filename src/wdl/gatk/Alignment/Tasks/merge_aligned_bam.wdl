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
   	   source ${BashPreamble}
   	   /bin/bash ${MergeBamScript} -b ${sep=',' InputBams} -s ${SampleName} -S ${SamtoolsExe} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      s_vmem: "${MergeSoftMemLimit}"
      h_vmem: "${MergeHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }
}
