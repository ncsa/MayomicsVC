#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/DesignBlock_2a/Tasks/realignment.wdl" as REALIGNMENT
import "src/wdl_scripts/DesignBlock_2a/Tasks/bqsr.wdl" as BQSR
import "src/wdl_scripts/DesignBlock_2a/Tasks/haplotyper.wdl" as HAPLOTYPER

workflow CallBlock2aTasks {

############## BOILERPLATE FOR DESIGN BLOCK 2a ########################################

   File GlobalAlignedSortedDedupedBam
   File GlobalAlignedSortedDedupedBamIdx

   String SampleName                                   

   File RefFasta                            
   File RefFai                                         

   File KnownSites                                     
   File KnownSitesIdx                                  
  
   File DBSNP
   File DBSNPIdx
   
   String Threads                                      
   String SentieonLicense                              
   String Sentieon                                     
   Boolean DebugMode                                   

   File RealignmentScript                              
   File BqsrScript
   File HaplotyperScript

######################################################################################
   
   call REALIGNMENT.realignmentTask  as realign {
      input:
         InputAlignedSortedDedupedBam = GlobalAlignedSortedDedupedBam,
         InputAlignedSortedDedupedBamIdx = GlobalAlignedSortedDedupedBamIdx,
         SampleName = SampleName,
         RefFasta = RefFasta,
         RefFai = RefFai,
         KnownSites = KnownSites,
         KnownSitesIdx = KnownSitesIdx,
         Threads = Threads,
         SentieonLicense = SentieonLicense,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         RealignmentScript = RealignmentScript
   }
   
   call BQSR.bqsrTask as bqsr {
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamIdx = realign.AlignedSortedDedupedRealignedBamIdx,
         SampleName = SampleName,
         RefFasta = RefFasta,
         RefFai = RefFai,
         KnownSites = KnownSites,
         KnownSitesIdx = KnownSitesIdx,
         DBSNP = DBSNP,
         DBSNPIdx = DBSNPIdx,
         Threads = Threads,
         SentieonLicense = SentieonLicense,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         BqsrScript = BqsrScript
   }

   call HAPLOTYPER.variantCallingTask as haplotype { 
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBamIdx,
         RecalTable = bqsr.RecalTable,
         SampleName = SampleName,
         RefFasta = RefFasta,
         DBSNP = DBSNP,
         DBSNPIdx = DBSNPIdx,
         Threads = Threads,
         SentieonLicense = SentieonLicense,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         HaplotyperScript = HaplotyperScript
   }
   
}
