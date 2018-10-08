##############################################################

###       This WDL script calls the CutAdapt WDL Task       ##

##############################################################

import "src/wdl_scripts/DesignBlock_1/Tasks/trim_sequences.wdl" as TRIMSEQ

workflow RunTrimInputSequencesTask {

   # tab-separated values, one lane per line
   # InputRead1, InputRead2
   #File InputReadsList
   #Array[Array[String]] InputReads = read_tsv(InputReadsList)
   Array[String]+ InputRead1
   Array[String]? InputRead2

   Int NumOfInputs = length(InputRead1)

   File Adapters                   # Adapter FastA File         
 
   String CutAdapt                 # Path to CutAdapt Tool
   String CutAdaptThreads          # Specifies the number of thread required per run

   Boolean PairedEnd               # Variable to check if single ended or not
   String DebugMode                # Variable to check if Debud Mode is on or not

   File TrimSeqScript              # Bash script which is called inside the WDL script
   File TrimEnvProfile             # File containing the environmental profile variables

   String SampleName               # Name of the Sample

   scatter (idx in range(NumOfInputs)) {
      # If PairedEnd=False, set InputRead2="null"
      if(PairedEnd) {
      call TRIMSEQ.trimsequencesTask as TRIMSEQ_paired {
         input:
            SampleName=SampleName,
            InputRead1=InputRead1[idx],
            InputRead2=InputRead2[idx],
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
            InputRead1=InputRead1[idx],
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
