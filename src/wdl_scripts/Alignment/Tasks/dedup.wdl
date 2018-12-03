###########################################################################################
##              This WDL script marks the duplicates on input sorted BAMs                ##
###########################################################################################

task dedupTask {

   File InputBams                  # Input Sorted BAM File
   File InputBais                  # Input Sorted Bam Index File

   String SampleName               # Name of the Sample

   String Sentieon                 # Variable path to Sentieon 

   String SentieonThreads          # Specifies the number of thread required per run
   String DebugMode                # Variable to check whether Debud Mode is on

   String DedupSoftMemLimit        # Soft memory limit - nice shutdown
   String DedupHardMemLimit        # Hard memory limit - kill immediately
   File BashPreamble               # shell file to source before each task
   File BashSharedFunctions        # Bash script with shared functions

   File DedupScript                # Bash script that is called inside the WDL script
   File DedupEnvProfile            # File containing the environmental profile variables

   command <<<
   	   source ${BashPreamble}
   	   /bin/bash ${DedupScript} -b ${InputBams} -s ${SampleName} -S ${Sentieon} -t ${SentieonThreads} -e ${DedupEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${DedupSoftMemLimit}"
      h_vmem: "${DedupHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.aligned.sorted.deduped.bam"
      File OutputBais = "${SampleName}.aligned.sorted.deduped.bam.bai"
   }
}
