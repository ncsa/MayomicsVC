###############################################################################################
####              This WDL script is used to run Alignment and HaplotyperVC blocks together  ##
###############################################################################################

import "src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/Alignment/TestTasks/Runalignment.wdl" as ALIGNMENT
import "src/wdl_scripts/Alignment/Tasks/dedup.wdl" as DEDUP

import "src/wdl_scripts/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment

import "src/wdl_scripts/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "src/wdl_scripts/HaplotyperVC/Tasks/bqsr.wdl" as BQSR
import "src/wdl_scripts/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER
import "src/wdl_scripts/HaplotyperVC/Tasks/vqsr.wdl" as VQSR

import "src/wdl_scripts/DeliveryOfHaplotyperVC/Tasks/deliver_HaplotyperVC.wdl" as DELIVER_HaplotyperVC


workflow MasterWF {

   call CUTADAPTTRIM.RunTrimSequencesTask as trimseq

   call ALIGNMENT.RunAlignmentTask as align {
      input:
         InputReads = trimseq.TrimmedInputReads
   }

   call DEDUP.dedupTask as dedup {
      input:
         InputAlignedSortedBam  = align.AlignedSortedBams,
         InputAlignedSortedBamBai = align.AlignedSortedBamBais
   }



   call DELIVER_Alignment.deliverAlignmentTask as DAB {
      input:
         AlignedSortedDedupedBam = dedup.AlignedSortedDedupedBam,
         AlignedSortedDedupedBamBai = dedup.AlignedSortedDedupedBamBai
   }



   call REALIGNMENT.realignmentTask  as realign {
      input:
         InputAlignedSortedDedupedBam = dedup.AlignedSortedDedupedBam,
         InputAlignedSortedDedupedBamBai = dedup.AlignedSortedDedupedBamBai
   }


   call BQSR.bqsrTask as bqsr {
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamBai = realign.AlignedSortedDedupedRealignedBamBai,
   }

   call HAPLOTYPER.variantCallingTask as haplotype {
      input:
         InputAlignedSortedDedupedRealignedBam = realign.AlignedSortedDedupedRealignedBam,
         InputAlignedSortedDedupedRealignedBamBai = realign.AlignedSortedDedupedRealignedBamBai,
         RecalTable = bqsr.RecalTable,
   }

   call VQSR.vqsrTask as vqsr {
      input:
         InputVCF = haplotype.VCF,
         InputVCFIdx = haplotype.VcfIdx,
   }



   call DELIVER_HaplotyperVC.deliverHaplotyperVCTask as DHVC {
      input:
         RecalibratedVcf = vqsr.RecalibratedVcf,
         RecalibratedVcfIdx = vqsr.RecalibratedVcfIdx
   } 

}






##################################################################################################

#import "src/wdl_scripts/Alignment/Workflow/Alignment_ForWorkflow.wdl" as WF1
#import "src/wdl_scripts/DeliveryOfAlignment/Workflow/DeliverAlignment_ForWorkflow.wdl" as DWF1
#import "src/wdl_scripts/HaplotyperVC/Workflow/HaplotyperVC_ForWorkflow.wdl" as WF2a
#import "src/wdl_scripts/DeliveryOfHaplotyperVC/Workflow/DeliverHaplotyperVC_ForWorkflow.wdl" as DWF2a
#
#workflow MasterWF {
#
#   call WF1.CallAlignmentTasks as AlignmentBlock
#
#   call DWF1.CallDeliveryAlignmentTask as DAB {
#      input:
#         GlobalAlignedSortedDedupedBam = AlignmentBlock.GlobalAlignedSortedDedupedBam,
#         GlobalAlignedSortedDedupedBamBai = AlignmentBlock.GlobalAlignedSortedDedupedBamBai
#   }
#
#   call WF2a.CallHaplotyperVCTasks as HaplotyperVCBlock {
#      input:
#         GlobalAlignedSortedDedupedBam = AlignmentBlock.GlobalAlignedSortedDedupedBam,
#         GlobalAlignedSortedDedupedBamBai = AlignmentBlock.GlobalAlignedSortedDedupedBamBai
#   }
#
#   call DWF2a.CallDeliveryHaplotyperVCTask as DHVCB {
#      input:
#         GlobalRecalibratedVcf = HaplotyperVCBlock.GlobalRecalibratedVcf,
#         GlobalRecalibratedVcfIdx = HaplotyperVCBlock.GlobalRecalibratedVcfIdx
#   }
#}
