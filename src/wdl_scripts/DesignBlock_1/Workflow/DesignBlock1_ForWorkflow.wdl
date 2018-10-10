##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "../TestTasks/Testtrim_sequences.wdl" as CUTADAPTTRIM
import "../TestTasks/Testalignment.wdl" as ALIGNMENT
import "../Tasks/dedup.wdl" as DEDUP 

workflow CallBlock1Tasks {

############## BOILERPLATE FOR DESIGN BLOCK 1 #######################################

   Array[Array[String]] InputReads
   File Adapters  
   String CutAdapt 
   String CutAdaptThreads    
   Boolean PairedEnd               
   String DebugMode
   String SampleName

   File Ref
   File RefAmb
   File RefAnn
   File RefBwt
   File RefPac
   File RefSa  
   String Group
   String Platform

   String Sentieon 
   String SentieonThreads    
   String ChunkSizeInBases

   File TrimSeqScript
   File TrimEnvProfile

   File AlignmentScript
   File AlignEnvProfile

   File DedupScript
   File DedupEnvProfile


#####################################################################################          
   
   call CUTADAPTTRIM.RunTrimInputSequencesTask as trimseq {
      input:
         InputReads = InputReads,
         Adapters = Adapters,
         CutAdapt = CutAdapt,
         CutAdaptThreads = CutAdaptThreads,
         PairedEnd = PairedEnd,
         DebugMode = DebugMode,
         TrimSeqScript = TrimSeqScript,
         TrimEnvProfile = TrimEnvProfile,
         SampleName = SampleName
   }
    
   call ALIGNMENT.CallalignmentTask as align {
      input:
         InputReads = trimseq.TrimmedInputReads,
         Ref = Ref,
         SampleName = SampleName,
         RefAmb = RefAmb,
         RefAnn = RefAnn,
         RefBwt = RefBwt,
         RefPac = RefPac,
         RefSa = RefSa,
         Sentieon = Sentieon,
         ChunkSizeInBases = ChunkSizeInBases,
         Group = Group,
         Platform = Platform,
         DebugMode = DebugMode,
         SentieonThreads = SentieonThreads,
         PairedEnd = PairedEnd,
         AlignmentScript = AlignmentScript,
         AlignEnvProfile = AlignEnvProfile
   }
   
   call DEDUP.dedupTask as dedup {
      input:
         InputAlignedSortedBam  = align.AlignedSortedBams,
         InputAlignedSortedBamBai = align.AlignedSortedBamBais,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         SentieonThreads = SentieonThreads,
         SampleName = SampleName,
         DedupScript = DedupScript,
         DedupEnvProfile = DedupEnvProfile
   }
    
   output {
     
      File GlobalAlignedSortedDedupedBam = dedup.AlignedSortedDeduppedBam
      File GlobalAlignedSortedDedupedBamBai = dedup.AlignedSortedDeduppedBamBai
   }    

}
