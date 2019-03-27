###############################################################################################
####              This WDL script is used to run Alignment and HaplotyperVC blocks together  ##
###############################################################################################

import "MayomicsVC/src/wdl/gatk/Alignment/TestTasks/Runtrim_sequences.wdl" as CUTADAPTTRIM
import "MayomicsVC/src/wdl/gatk/Alignment/TestTasks/Runalignment.wdl" as ALIGNMENT
import "MayomicsVC/src/wdl/gatk/Alignment/Tasks/merge_aligned_bam.wdl" as MERGEBAM
import "MayomicsVC/src/wdl/gatk/Alignment/Tasks/dedup.wdl" as DEDUP

import "MayomicsVC/src/wdl/gatk/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment

import "MayomicsVC/src/wdl/gatk/HaplotyperVC/Tasks/bqsr.wdl" as BQSR
import "MayomicsVC/src/wdl/gatk/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER

import "MayomicsVC/src/wdl/gatk/HaplotyperVC/Tasks/merge_gvcfs.wdl" as MERGEGVCF
import "MayomicsVC/src/wdl/gatk/DeliveryOfHaplotyperVC/Tasks/deliver_HaplotyperVC.wdl" as DELIVER_HaplotyperVC

workflow GermlineMasterWF {

   Array[Array[File]] NormalInputReads

   Boolean Trimming
   Boolean MarkDuplicates
   
   Array[String] GenomicIntervals

   if(Trimming) {

      call CUTADAPTTRIM.RunTrimSequencesTask as trimseq {
         input:
            InputReads = NormalInputReads
      }
   }
   
   Array[Array[File]] AlignInputReads = select_first([trimseq.Outputs,NormalInputReads])


   call ALIGNMENT.RunAlignmentTask as align {
      input:
         InputReads = AlignInputReads
   }

   call MERGEBAM.mergebamTask as merge {
      input:
         InputBams = align.OutputBams,
         InputBais = align.OutputBais
   }

   if(MarkDuplicates) {
   
      call DEDUP.dedupTask as dedup {
         input:
            InputBams = merge.OutputBams,
            InputBais = merge.OutputBais
      }
   }

   File DeliverAlignOutputBams = select_first([dedup.OutputBams,merge.OutputBams])
   File DeliverAlignOutputBais = select_first([dedup.OutputBais,merge.OutputBais])


   call DELIVER_Alignment.deliverAlignmentTask as DAB {
      input:
         InputBams = DeliverAlignOutputBams,
         InputBais = DeliverAlignOutputBais
   }



   scatter (interval in GenomicIntervals) {

      call BQSR.bqsrTask as bqsr {
         input:
            InputBams = DeliverAlignOutputBams,
            InputBais = DeliverAlignOutputBais,
            GenomicInterval = interval
      }
   
      call HAPLOTYPER.variantCallingTask as haplotype {
         input:
            InputBams = bqsr.OutputBams,
            InputBais = bqsr.OutputBais,
            GenomicInterval = interval
      }

   }

   call MERGEGVCF.mergegvcfsTask as mergegvcfs {
      input:
         InputGvcfs = haplotype.OutputVcf,
         InputIdxs = haplotype.OutputVcfIdx
   }
  
   call DELIVER_HaplotyperVC.deliverHaplotyperVCTask as DHVC {
      input:
         InputVcf = mergegvcfs.OutputVcf,
         InputVcfIdx = mergegvcfs.OutputVcfIdx
    } 


}
