##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/Alignment/Tasks/trim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/Alignment/Tasks/alignment.wdl" as ALIGNMENT
import "src/wdl_scripts/Alignment/Tasks/dedup.wdl" as DEDUP 

workflow CallAlignmentTasks {

############## BOILERPLATE FOR DESIGN BLOCK 1 #######################################

   File InputRead1 
   String InputRead2
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
   String BWAExtraOptionsString

   File TrimSeqScript
   File TrimEnvProfile

   File AlignmentScript
   File AlignEnvProfile

   File DedupScript
   File DedupEnvProfile


#####################################################################################          
   
   call CUTADAPTTRIM.trimsequencesTask as trimseq {
      input:
         InputRead1 = InputRead1,
         InputRead2 = InputRead2,
         Adapters = Adapters,
         CutAdapt = CutAdapt,
         CutAdaptThreads = CutAdaptThreads,
         PairedEnd = PairedEnd,
         DebugMode = DebugMode,
         TrimSeqScript = TrimSeqScript,
         TrimEnvProfile = TrimEnvProfile,
         SampleName = SampleName
   }
    
   call ALIGNMENT.alignmentTask as align {
      input:
         InputRead1 = trimseq.TrimmedInputRead1,
         InputRead2 = trimseq.TrimmedInputRead2,
         Ref = Ref,
         SampleName = SampleName,
         RefAmb = RefAmb,
         RefAnn = RefAnn,
         RefBwt = RefBwt,
         RefPac = RefPac,
         RefSa = RefSa,
         Sentieon = Sentieon,
         ChunkSizeInBases = ChunkSizeInBases,
         BWAExtraOptionsString = BWAExtraOptionsString,
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
         InputAlignedSortedBam  = align.AlignedSortedBam,
         InputAlignedSortedBamBai = align.AlignedSortedBamBai,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         SentieonThreads = SentieonThreads,
         SampleName = SampleName,
         DedupScript = DedupScript,
         DedupEnvProfile = DedupEnvProfile
   }
    
   output {
     
      File GlobalAlignedSortedDedupedBam = dedup.AlignedSortedDedupedBam
      File GlobalAlignedSortedDedupedBamBai = dedup.AlignedSortedDedupedBamBai
   }    

}
