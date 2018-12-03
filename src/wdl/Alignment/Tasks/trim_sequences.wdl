###########################################################################################

##              This WDL scripts trim the Inputs Fasta File using CutAdapt               ##

###########################################################################################         

task trimsequencesTask {

   File InputRead1                 # Input Read File             
   String InputRead2               # Input Read File             

   String SampleName               # Name of the Sample

   File Adapters                   # Adapter FastA File         
 
   String CutAdapt                 # Path to CutAdapt Tool
   String CutAdaptThreads          # Number of threads for cutadapt to use

   Boolean PairedEnd               # Variable to check if single ended or not

   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File TrimSeqScript              # Bash script which is called inside the WDL script
   File TrimEnvProfile             # File containing the environmental profile variables

   String TrimSoftMemLimit         # Soft memory limit - nice shutdown
   String TrimHardMemLimit         # Hard memory limit - kill immediately

   String DebugMode                # Variable to check if Debug Mode is on or not


   command <<<
     source ${BashPreamble}
     /bin/bash ${TrimSeqScript} -P ${PairedEnd} -l ${InputRead1} -r ${InputRead2} -s ${SampleName} -A ${Adapters} -C ${CutAdapt} -t ${CutAdaptThreads} -e ${TrimEnvProfile} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${CutAdaptThreads}"
      s_vmem: "${TrimSoftMemLimit}"
      h_vmem: "${TrimHardMemLimit}"
   }

   output {
      Array[File] Outputs = glob("${SampleName}.read?.trimmed.fq.gz")
   }
}

