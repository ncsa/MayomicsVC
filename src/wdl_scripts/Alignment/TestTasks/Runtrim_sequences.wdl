##############################################################

###       This WDL script calls the CutAdapt WDL Task       ##

##############################################################

import "src/wdl_scripts/Alignment/Tasks/trim_sequences.wdl" as TRIMSEQ

workflow RunTrimSequencesTask {

   Array[Array[File]] InputReads   # One lane per subarray with one or two input reads

   File Adapters                   # Adapter FastA File         
 
   String CutAdapt                 # Path to CutAdapt Tool
   String CutAdaptThreads          # Specifies the number of thread required per run

   Boolean PairedEnd               # Variable to check if single ended or not
   String DebugMode                # Variable to check if Debud Mode is on or not

   File TrimSeqScript              # Bash script which is called inside the WDL script
   File TrimEnvProfile             # File containing the environmental profile variables

   String SampleName               # Name of the Sample


   scatter (lane in InputReads) {
      # If PairedEnd=False, set InputRead2="null"
      if(PairedEnd) {
         call TRIMSEQ.trimsequencesTask as TRIMSEQ_paired {
            input:
               SampleName=SampleName,
               InputRead1=lane[0],
               InputRead2=lane[1],
               Adapters=Adapters,
               CutAdapt=CutAdapt,
               CutAdaptThreads=CutAdaptThreads,
               PairedEnd=PairedEnd,
               DebugMode=DebugMode,
               TrimSeqScript=TrimSeqScript,
               TrimEnvProfile=TrimEnvProfile
         }
      }

      if(!PairedEnd) {
         call TRIMSEQ.trimsequencesTask as TRIMSEQ_single {
            input:
               SampleName=SampleName,
               InputRead1=lane[0],
               InputRead2="null",
               Adapters=Adapters,
               CutAdapt=CutAdapt,
               CutAdaptThreads=CutAdaptThreads,
               PairedEnd=PairedEnd,
               DebugMode=DebugMode,
               TrimSeqScript=TrimSeqScript,
               TrimEnvProfile=TrimEnvProfile
         }
      }
   }

   output {
      # Unify outputs from scatter and filter out null entries 
      Array[Array[File]] Outputs = select_all(flatten([TRIMSEQ_paired.Outputs,TRIMSEQ_single.Outputs]))
   }
} 
