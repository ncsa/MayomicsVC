###########################################################################################

##              This WDL script performs realignment using Sentieon                      ##

###########################################################################################

task realignmentTask {

   File InputBams                                     # Input Sorted Deduped Bam
   File InputBais                                     # Input Sorted Deduped Bam Index

   File Ref                                           # Reference Genome
   File RefFai                                        # Reference Index File
                                  
   String SampleName                                  # Name of the Sample

   String RealignmentKnownSites                       # List of known sites

   String Sentieon                                    # Path to Sentieon
   String SentieonThreads                             # No of Threads for the Tool

   String RealignSoftMemLimit                         # Soft memory limit - nice shutdown
   String RealignHardMemLimit                         # Hard memory limit - kill immediately

   String DebugMode                                   # Enable or Disable Debug Mode

   File BashPreamble                                  # Bash script that helps control zombie processes
   File BashSharedFunctions                           # Bash script that contains shared helpful functions
   File RealignmentScript                             # Path to bash script called within WDL script
   File RealignEnvProfile                             # File containing the environmental profile variables

   command <<<
      source ${BashPreamble}
      /bin/bash ${RealignmentScript} -s ${SampleName} -b ${InputBams} -G ${Ref} -k ${RealignmentKnownSites} -S ${Sentieon} -t ${SentieonThreads} -e ${RealignEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${RealignSoftMemLimit}"
      h_vmem: "${RealignHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.aligned.sorted.deduped.realigned.bam"
      File OutputBais = "${SampleName}.aligned.sorted.deduped.realigned.bam.bai"
   }  
} 
