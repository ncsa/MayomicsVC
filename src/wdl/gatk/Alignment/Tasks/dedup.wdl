###########################################################################################
##              This WDL script marks the duplicates on input sorted BAMs                ##
###########################################################################################

task dedupTask {
   
   String SampleName               # Name of the Sample

   File InputBams                  # Input Sorted BAM File
   File InputBais                  # Input Sorted Bam Index File

   
   File GATKExe                    # Path to GATK4 executable
   File JavaExe                    # Path to Java8 executable
   String JavaOptionsString        # String of java vm options, like garbage collection and maximum and minimum memory. Can NOT be empty

   String DebugMode                # Variable to check whether Debud Mode is on

   String DedupSoftMemLimit        # Soft memory limit - nice shutdown
   String DedupHardMemLimit        # Hard memory limit - kill immediately
   File BashPreamble               # shell file to source before each task
   File BashSharedFunctions        # Bash script with shared functions

   File DedupScript                # Bash script that is called inside the WDL script

   command <<<
   	   source ${BashPreamble}
   	   /bin/bash ${DedupScript} -s ${SampleName} -b ${InputBams} -S ${GATKExe} -J ${JavaExe} -e ${JavaOptionsString} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: 1
      s_vmem: "${DedupSoftMemLimit}"
      memory: "${DedupHardMemLimit}"
      docker : "broadinstitute/gatk:4.1.1.0"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }
}
