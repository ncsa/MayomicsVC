#################################################################################################
#                                                                                               #  
##              This WDL script merges the input asligned BAMs before deduplication            ##
#                                                                                               #  
#################################################################################################

task mergebamTask {

   Array[File] InputBams           # Input Sorted BAM File
   Array[File] InputBais           # Input Sorted Bam Index File

   String SampleName               # Name of the Sample

   String Sentieon                 # Variable path to Sentieon 

   String SentieonThreads          # Specifies the number of thread required per run
   String DebugMode                # Variable to check whether Debud Mode is on

   String MergeSoftMemLimit        # Soft memory limit - nice shutdown
   String MergeHardMemLimit        # Hard memory limit - kill immediately
   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File MergeBamScript             # Bash script that is called inside the WDL script
   File MergeBamEnvProfile         # File containing the environmental profile variables

   command <<<
   	   source ${BashPreamble}
   	   /bin/bash ${MergeBamScript} -b ${sep=',' InputBams} -s ${SampleName} -S ${Sentieon} -t ${SentieonThreads} -e ${MergeBamEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${MergeSoftMemLimit}"
      h_vmem: "${MergeHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.aligned.sorted.merged.bam"
      File OutputBais = "${SampleName}.aligned.sorted.merged.bam.bai"
   }
}
