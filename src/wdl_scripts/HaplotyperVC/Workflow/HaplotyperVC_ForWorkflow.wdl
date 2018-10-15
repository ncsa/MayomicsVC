#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "src/wdl_scripts/HaplotyperVC/Tasks/bqsr.wdl" as BQSR
import "src/wdl_scripts/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER
import "src/wdl_scripts/HaplotyperVC/Tasks/vqsr.wdl" as VQSR

workflow CallHaplotyperVCTasks {

   File GlobalAlignedSortedDedupedBam
   File GlobalAlignedSortedDedupedBamBai

############## BOILERPLATE FOR DESIGN BLOCK 2a ########################################

   String SampleName                                   

   File Ref                            
   File RefFai                                         

   String RealignmentKnownSites                                     
   String BQSRKnownSites                                     
   File DBSNP
   File DBSNPIdx

   String VqsrSnpResourceString
   String VqsrIndelResourceString
   String AnnotateText
  
   String Sentieon                                     
   String SentieonThreads                                      

   String DebugMode                                   

   File RealignmentScript                              
   File RealignEnvProfile

   File BqsrScript
   File BqsrEnvProfile

   File HaplotyperScript
   File HaplotyperEnvProfile

   String HaplotyperExtraOptions

   File VqsrScript
   File VqsrEnvProfile

######################################################################################
   
   call REALIGNMENT.realignmentTask  as realign {
      input:
         InputAlignedSortedDedupedBam = GlobalAlignedSortedDedupedBam,
         InputAlignedSortedDedupedBamBai = GlobalAlignedSortedDedupedBamBai,
         SampleName = SampleName,
         Ref = Ref,
         RefFai = RefFai,
         RealignmentKnownSites = RealignmentKnownSites,
         SentieonThreads = SentieonThreads,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         RealignmentScript = RealignmentScript,
         RealignEnvProfile = RealignEnvProfile
   }
   
   call BQSR.bqsrTask as bqsr {
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamBai = realign.AlignedSortedDedupedRealignedBamBai,
         SampleName = SampleName,
         Ref = Ref,
         RefFai = RefFai,
         BQSRKnownSites = BQSRKnownSites,
         SentieonThreads = SentieonThreads,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         BqsrScript = BqsrScript,
         BqsrEnvProfile = BqsrEnvProfile
   }


   call HAPLOTYPER.variantCallingTask as haplotype { 
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamBai = realign.AlignedSortedDedupedRealignedBamBai,
         RecalTable = bqsr.RecalTable,
         SampleName = SampleName,
         Ref = Ref,
         RefFai = RefFai,
         DBSNP = DBSNP,
         DBSNPIdx = DBSNPIdx,
         SentieonThreads = SentieonThreads,
         Sentieon = Sentieon,
         DebugMode = DebugMode,
         HaplotyperScript = HaplotyperScript,
         HaplotyperEnvProfile = HaplotyperEnvProfile,
	 HaplotyperExtraOptions = HaplotyperExtraOptions
   }

   call VQSR.vqsrTask as vqsr {
      input:
         InputVCF = haplotype.VCF,
         InputVCFIdx = haplotype.VcfIdx,
         SampleName = SampleName,
         Ref = Ref,
         RefFai = RefFai,
         VqsrSnpResourceString=VqsrSnpResourceString,
         VqsrIndelResourceString=VqsrIndelResourceString,
         SentieonThreads = SentieonThreads,
         Sentieon = Sentieon,
         AnnotateText = AnnotateText,
         DebugMode = DebugMode,
         VqsrScript = VqsrScript,
         VqsrEnvProfile = VqsrEnvProfile
   }

   output {
      File GlobalRecalibratedVcf = vqsr.RecalibratedVcf
      File GlobalRecalibratedVcfIdx = vqsr.RecalibratedVcfIdx
   }

}
