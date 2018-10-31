###########################################################################################

##              This WDL scripts trim the Inputs Fasta File using CutAdapt               ##

##                                    Script Options                     
#         -t        "Number of Threads"                         (Required)
#         -P        "Single Ended Reads specification"          (Required)
#         -r        "Left Fastq File"                           (Required)
#         -R        "Right Fastq File"                          (Optional)
#         -s        "Name of the sample"                        (Optional)
#         -A        "Adapter File for CutAdapt"                 (Required)
#         -C        "Path to CutAdapt Tool"                     (Required)
#         -e        "Path to the environmental profile          (Required)
#         -d        "debug mode on/off                          (Optional: can be empty)

###########################################################################################         

task trimsequencesTask {

   File InputRead1                 # Input Read File             
   String InputRead2               # Input Read File             

   String SampleName               # Name of the Sample

   File Adapters                   # Adapter FastA File         
 
   String CutAdapt                 # Path to CutAdapt Tool
   String CutAdaptThreads          # Number of threads for cutadapt to use

   Boolean PairedEnd               # Variable to check if single ended or not

   File TrimSeqScript              # Bash script which is called inside the WDL script
   File TrimEnvProfile             # File containing the environmental profile variables

   String SoftMemLimit             # Soft memory limit - nice shutdown
   String HardMemLimit             # Hard memory limit - kill immediately

   String DebugMode                # Variable to check if Debug Mode is on or not


   command {
      /bin/bash ${TrimSeqScript} -P ${PairedEnd} -l ${InputRead1} -r ${InputRead2} -s ${SampleName} -A ${Adapters} -C ${CutAdapt} -t ${CutAdaptThreads} -e ${TrimEnvProfile} ${DebugMode}
   }

   runtime {
      cpu: "${CutAdaptThreads}"
      s_vmem: "${SoftMemLimit}"
      h_vmem: "${HardMemLimit}"
   }

   output {
      Array[File] Outputs = glob("${SampleName}.read?.trimmed.fq.gz")
   }

}

