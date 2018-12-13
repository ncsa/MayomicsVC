###############################################################################################
####                  Standalone script for realignment and variant calling                  ##
###############################################################################################

import "MayomicsVC/src/wdl/HaplotyperVC/Tasks/realignment.wdl" as REALIGNMENT
import "MayomicsVC/src/wdl/SomaticVC/Tasks/strelka.wdl" as STRELKA
import "MayomicsVC/src/wdl/SomaticVC/Tasks/mutect.wdl" as MUTECT
import "MayomicsVC/src/wdl/SomaticVC/Tasks/combine_variants.wdl" as MERGEVCF

import "MayomicsVC/src/wdl/DeliveryOfSomaticVC/Tasks/deliver_SomaticVC.wdl" as DELIVER_SomaticVC


workflow SomaticMasterWF {

   File NormalInputBams
   File NormalInputBais
   File TumorInputBams
   File TumorInputBais
   String SampleNameNormal
   String SampleNameTumor

   call REALIGNMENT.realignmentTask  as TumorRealign {
      input:
         SampleName = SampleNameTumor,
         InputBams = TumorInputBams,
         InputBais = TumorInputBais
   }

   call REALIGNMENT.realignmentTask  as NormalRealign {
      input:
         SampleName = SampleNameNormal,
         InputBams = NormalInputBams,
         InputBais = NormalInputBais
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
         InputVcf = merge_somatic_vcf.OutputVcf,
         InputVcfIdx = merge_somatic_vcf.OutputVcfIdx
   } 

}
