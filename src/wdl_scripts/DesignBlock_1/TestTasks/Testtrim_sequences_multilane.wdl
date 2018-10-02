##############################################################

###       This WDL script calls the CutAdapt WDL Task       ##

##############################################################

import "src/wdl_scripts/DesignBlock_1/Tasks/trim_sequences.wdl" as TRIMSEQ

workflow RunTrimInputSequencesTask {

   # tab-separated values, one lane per line
   # InputRead1, InputRead2
   #File InputReadsList
   #Array[Array[String]] InputReads = read_tsv(InputReadsList)
   Array[String] InputRead1
   Array[String] InputRead2

   Array[Pair[String,String]] InputReads = zip(InputRead1, InputRead2)

   File Adapters                   # Adapter FastA File         
 
   String CutAdapt                 # Path to CutAdapt Tool
   String CutAdaptThreads          # Specifies the number of thread required per run

   Boolean PairedEnd               # Variable to check if single ended or not
   String DebugMode                # Variable to check if Debud Mode is on or not

   File TrimSeqScript              # Bash script which is called inside the WDL script
   File TrimEnvProfile             # File containing the environmental profile variables

   String SampleName               # Name of the Sample

   scatter (lane in InputReads) {
      # TODO: If PairedEnd=False, set InputRead2=NULL; otherwise verify InputRead2 is valid
      if(PairedEnd) {
      call TRIMSEQ.trimsequencesTask as TRIMSEQ_paired {
         input:
            SampleName=SampleName,
            InputRead1=lane.left,
            InputRead2=lane.right,
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
            InputRead1=lane.left,
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

} 
