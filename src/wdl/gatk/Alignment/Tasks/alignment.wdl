###########################################################################################
##              This WDL script performs alignment using BWA Mem                         ##
###########################################################################################

task alignmentTask {
   String SampleName               # Name of the Sample
   String Platform                 # sequencing platform for read group
   String Library                  # Sequencing library for read group
   String PlatformUnit             # Platform unit / flowcell ID for read group
   String CenterName               # Name of the sequencing center for read group
   Boolean PairedEnd               # Variable to check if single ended or not
   File InputRead1                 # Input Read File           
   String InputRead2               # Input Read File           
   File Ref                        # Reference Genome
   File RefAmb                     # reference file index
   File RefAnn                     # reference file index
   File RefBwt                     # reference file index
   File RefPac                     # reference file index
   File RefSa                      # reference file index
   File BWAExe                     # Path to BWA executable**************
   String ChunkSizeInBases         # The -K option for BWA MEM
   String BWAExtraOptionsString    # String of extra options for BWA. This can be an empty string.
   File SamtoolsExe                # Path to samtools executable
   String Threads                  # Specifies the number of thread required per run
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   String DebugMode                # Flag to enable Debug Mode

   File BashPreamble               # Bash script that helps control zombie processes
   File AlignmentScript            # Bash script which is called inside the WDL script
  
   String AlignSoftMemLimit        # Soft memory limit - nice shutdown
   String AlignHardMemLimit        # Hard memory limit - kill immediately


   command <<<
      source ${BashPreamble}
      /bin/bash ${AlignmentScript} -s ${SampleName} -p ${Platform} -L ${Library} -f ${PlatformUnit} -c ${CenterName} -P ${PairedEnd} -l ${InputRead1} -r ${InputRead2} -G ${Ref} -e ${BWAExe} -K ${ChunkSizeInBases} -o "'${BWAExtraOptionsString}'" -S ${SamtoolsExe} -t ${Threads} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${Threads}"
      s_vmem: "${AlignSoftMemLimit}"
      h_vmem: "${AlignHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }

} 

