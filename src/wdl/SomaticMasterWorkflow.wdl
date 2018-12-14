###############################################################################################
####              This WDL script is used to run Alignment and HaplotyperVC blocks together  ##
###############################################################################################

import "MayomicsVC/src/wdl/Alignment/TestTasks/Runtrim_sequences.wdl" as CUTADAPTTRIM
import "MayomicsVC/src/wdl/Alignment/TestTasks/Runalignment.wdl" as ALIGNMENT
import "MayomicsVC/src/wdl/Alignment/Tasks/merge_aligned_bam.wdl" as MERGEBAM
import "MayomicsVC/src/wdl/Alignment/Tasks/dedup.wdl" as DEDUP

import "MayomicsVC/src/wdl/DeliveryOfAlignment/Tasks/deliver_alignment.wdl" as DELIVER_Alignment


import "MayomicsVC/src/wdl/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "MayomicsVC/src/wdl/SomaticVC/Tasks/strelka.wdl" as STRELKA
import "MayomicsVC/src/wdl/SomaticVC/Tasks/mutect.wdl" as MUTECT
import "MayomicsVC/src/wdl/SomaticVC/Tasks/combine_variants.wdl" as MERGEVCF

import "MayomicsVC/src/wdl/DeliveryOfSomaticVC/Tasks/deliver_SomaticVC.wdl" as DELIVER_SomaticVC


workflow SomaticMasterWF {

   Array[Array[File]] NormalInputReads
   Array[Array[File]] TumorInputReads
   String SampleNameNormal
   String SampleNameTumor

   Boolean Trimming 


   ######     Alignment subworkflow for Tumor sample      ######
   #############################################################

   if(Trimming) {
      call CUTADAPTTRIM.RunTrimSequencesTask as TumorTrimseq {
         input: 
            InputReads = TumorInputReads 
      }
   }

   Array[Array[File]] TumorAlignInputReads = select_first([TumorTrimseq.Outputs,TumorInputReads])

   call ALIGNMENT.RunAlignmentTask as TumorAlign {
      input:
         SampleName = SampleNameTumor,
         InputReads = TumorAlignInputReads

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
         InputBais = TumorDedup.OutputBais,
   }




   ######     Alignment subworkflow for Normal sample      ######
   ##############################################################

   if(Trimming) {
      call CUTADAPTTRIM.RunTrimSequencesTask as NormalTrimseq {
         input:
            InputReads = NormalInputReads
      }
   }

   Array[Array[File]] NormalAlignInputReads = select_first([NormalTrimseq.Outputs,NormalInputReads])

   call ALIGNMENT.RunAlignmentTask as NormalAlign {
      input:
         SampleName = SampleNameNormal,
         InputReads = NormalAlignInputReads
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

   call MERGEVCF.combineVariantsTask as merge_somatic_vcf {
      input:
         StrelkaVcfBgz = strelka.OutputVcfBgz,
         StrelkaVcfBgzTbi = strelka.OutputVcfBgzTbi,
         MutectVcfBgz = mutect.OutputVcfBgz,
         MutectVcfBgzTbi = mutect.OutputVcfBgzTbi
   }

   call DELIVER_SomaticVC.deliverSomaticVCTask as DSVC {
      input:
         InputVcf = merge_somatic_vcf.OutputVcfGz,
         InputVcfIdx = merge_somatic_vcf.OutputVcfGzTbi
   } 

}
