###############################################################################################
####              This WDL script is used to run Alignment and HaplotyperVC blocks together  ##
###############################################################################################

import "src/wdl/Alignment/TestTasks/Runtrim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl/Alignment/TestTasks/Runalignment.wdl" as ALIGNMENT
import "src/wdl/Alignment/Tasks/merge_aligned_bam.wdl" as MERGEBAM
import "src/wdl/Alignment/Tasks/dedup.wdl" as DEDUP

import "src/wdl/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment


import "src/wdl/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "src/wdl/SomaticVC/Tasks/strelka.wdl" as STRELKA
import "src/wdl/SomaticVC/Tasks/mutect.wdl" as MUTECT
import "src/wdl/SomaticVC/Tasks/merge_somatic_vcf.wdl" as MERGEVCF

import "src/wdl/DeliveryOfSomaticVC/Tasks/deliver_SomaticVC.wdl" as DELIVER_SomaticVC


workflow SomaticMasterWF {

   Array[Array[File]] NormalInputReads
   Array[Array[File]] TumorInputReads
   String SampleNameNormal
   String SampleNameTumor


   ######     Alignment subworkflow for Tumor sample      ######
   #############################################################

   call CUTADAPTTRIM.RunTrimSequencesTask as TumorTrimseq {
      input: 
         SampleName = SampleNameTumor,
         InputReads = TumorInputReads 
   }


   call ALIGNMENT.RunAlignmentTask as TumorAlign {
      input:
         SampleName = SampleNameTumor,
         InputReads = TumorTrimseq.Outputs
   }

   call MERGEBAM.mergebamTask as TumorMergeBam {
      input:
         SampleName = SampleNameTumor,
         InputBams = TumorAlign.OutputBams,
         InputBais = TumorAlign.OutputBais
   }

   call DEDUP.dedupTask as TumorDedup {
      input:
         SampleName = SampleNameTumor,
         InputBams = TumorMergeBam.OutputBams,
         InputBais = TumorMergeBam.OutputBais
   }

   call DELIVER_Alignment.deliverAlignmentTask as TumorDAB {
      input:
         SampleName = SampleNameTumor,
         InputBams = TumorDedup.OutputBams,
         InputBais = TumorDedup.OutputBais
   }




   ######     Alignment subworkflow for Normal sample      ######
   ##############################################################

   call CUTADAPTTRIM.RunTrimSequencesTask as NormalTrimseq {
      input:
         SampleName = SampleNameNormal,
         InputReads = NormalInputReads
   }

   call ALIGNMENT.RunAlignmentTask as NormalAlign {
      input:
         SampleName = SampleNameNormal,
         InputReads = NormalTrimseq.Outputs
   }

   call MERGEBAM.mergebamTask as NormalMergeBam {
      input:
         SampleName = SampleNameNormal,
         InputBams = NormalAlign.OutputBams,
         InputBais = NormalAlign.OutputBais
   }

   call DEDUP.dedupTask as NormalDedup {
      input:
         SampleName = SampleNameNormal,
         InputBams = NormalMergeBam.OutputBams,
         InputBais = NormalMergeBam.OutputBais
   }

   call DELIVER_Alignment.deliverAlignmentTask as NormalDAB {
      input:
         SampleName = SampleNameNormal,
         InputBams = NormalDedup.OutputBams,
         InputBais = NormalDedup.OutputBais,
   }





   ######     Realignment and variant calling subworkflow`  ######
   ##############################################################
   

   call REALIGNMENT.realignmentTask  as TumorRealign {
      input:
         SampleName = SampleNameTumor,
         InputBams = TumorDedup.OutputBams,
         InputBais = TumorDedup.OutputBais
   }

   call REALIGNMENT.realignmentTask  as NormalRealign {
      input:
         SampleName = SampleNameNormal,
         InputBams = NormalDedup.OutputBams,
         InputBais = NormalDedup.OutputBais
   }

   call STRELKA.strelkaTask as strelka {
      input:
         TumorBams = TumorRealign.OutputBams,
         TumorBais = TumorRealign.OutputBais,
         NormalBams = NormalRealign.OutputBams,
         NormalBais = NormalRealign.OutputBais
   }

   call MUTECT.mutectTask as mutect {
      input:
         TumorBams = TumorRealign.OutputBams,
         TumorBais = TumorRealign.OutputBais,
         NormalBams = NormalRealign.OutputBams,
         NormalBais = NormalRealign.OutputBais
   }

   call MERGEVCF.mergeSomaticVcfTask as merge_somatic_vcf {
      input:
         InputStrelkaVcf = strelka.OutputVcf,
         InputStrelkaVcfIdx = strelka.OutputVcfIdx,
         InputMutectVcf = mutect.OutputVcf,
         InputMutectVcfIdx = mutect.OutputVcfIdx
   }

   call DELIVER_SomaticVC.deliverSomaticVCTask as DSVC {
      input:
         InputVcf = merge_somatic_vcf.OutputVcf,
         InputVcfIdx = merge_somatic_vcf.OutputVcfIdx,
   } 

}
