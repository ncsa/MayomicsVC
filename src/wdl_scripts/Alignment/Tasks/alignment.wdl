###########################################################################################
##              This WDL script performs alignment using BWA Mem                         ##
###########################################################################################

task alignmentTask {

   File InputRead1                 # Input Read File           
   String InputRead2               # Input Read File           
   String SampleName               # Name of the Sample
   String Platform                 # sequencing platform for read group
   String Library                  # Sequencing library for read group
   String PlatformUnit             # Platform unit / flowcell ID for read group
   String CenterName               # Name of the sequencing center for read group
   Boolean PairedEnd               # Variable to check if single ended or not

   File Ref                        # Reference Genome
   File RefAmb                     # reference file index
   File RefAnn                     # reference file index
   File RefBwt                     # reference file index
   File RefPac                     # reference file index
   File RefSa                      # reference file index

   String Sentieon                 # Path to Sentieon
   String SentieonThreads          # Specifies the number of thread required per run

   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File AlignmentScript            # Bash script which is called inside the WDL script
   File AlignEnvProfile            # File containing the environmental profile variables
   String ChunkSizeInBases         # The -K option for BWA MEM
   String BWAExtraOptionsString    # String of extra options for BWA. This can be an empty string.

   String AlignSoftMemLimit        # Soft memory limit - nice shutdown
   String AlignHardMemLimit        # Hard memory limit - kill immediately

   String DebugMode                # Flag to enable Debug Mode


   command <<<
      source ${BashPreamble}
      /bin/bash ${AlignmentScript} -P ${PairedEnd} -l ${InputRead1} -r ${InputRead2} -s ${SampleName} -p ${Platform} -L ${Library} -f ${PlatformUnit} -c ${CenterName} -G ${Ref} -o "'${BWAExtraOptionsString}'" -K ${ChunkSizeInBases} -S ${Sentieon} -t ${SentieonThreads} -e ${AlignEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${AlignSoftMemLimit}"
      h_vmem: "${AlignHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.aligned.sorted.bam"
      File OutputBais = "${SampleName}.aligned.sorted.bam.bai"
   }

} 

