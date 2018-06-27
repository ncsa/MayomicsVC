##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/DesignBlock_1/Tasks/trim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/DesignBlock_1/Tasks/alignment.wdl" as ALIGNMENT
import "src/wdl_scripts/DesignBlock_1/Tasks/dedup.wdl" as DEDUP 

workflow CallBlock1Tasks {

############## BOILERPLATE FOR DESIGN BLOCK 1 #######################################

   File InputRead1 
   String InputRead2
   File Adapters  
   String CutAdapt 
   String Threads    
   Boolean PairedEnd               
   String DebugMode
   String SampleName

   File RefFasta
   File RefAmb
   File RefAnn
   File RefBwt
   File RefPac
   File RefSa  
   String Group
   String Platform

   String SentieonLicense 
   String Sentieon 

   File TrimSeqScript
   File AlignmentScript
   File DedupScript

#####################################################################################          
   
   call CUTADAPTTRIM.trimsequencesTask as trimseq {
      input:
         InputRead1 = InputRead1,
         InputRead2 = InputRead2,
         Adapters = Adapters,
         CutAdapt = CutAdapt,
         Threads = Threads,
         PairedEnd = PairedEnd,
         DebugMode = DebugMode,
         TrimSeqScript = TrimSeqScript,
         SampleName = SampleName
   }
    
   call ALIGNMENT.alignmentTask as align {
      input:
         InputRead1 = trimseq.TrimmedInputRead1,
         InputRead2 = trimseq.TrimmedInputRead2,
         RefFasta = RefFasta,
         SampleName = SampleName,
         RefAmb = RefAmb,
         RefAnn = RefAnn,
         RefBwt = RefBwt,
         RefPac = RefPac,
         RefSa = RefSa,
         SentieonLicense = SentieonLicense,
         Sentieon = Sentieon,
         Group = Group,
         Platform = Platform,
         DebugMode = DebugMode,
         Threads = Threads,
         PairedEnd = PairedEnd,
         AlignmentScript = AlignmentScript
   }
   
   call DEDUP.dedupTask as dedup {
      input:
         InputAlignedSortedBam  = align.AlignedSortedBam,
         InputAlignedSortedBamIdx = align.AlignedSortedBamIdx,
         SentieonLicense = SentieonLicense,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         Threads = Threads,
         SampleName = SampleName,
         DedupScript = DedupScript
   }
    
   output {
     
      File GlobalAlignedSortedDedupedBam = dedup.AlignedSortedDeduppedBam
      File GlobalAlignedSortedDedupedBamIdx = dedup.AlignedSortedDeduppedBamIdx
   }    

}
