#########################################################################################################
####              This WDL script is used to run the  steps as individual modules              ##
##########################################################################################################

import "MayomicsVC/src/wdl/gatk/HaplotyperVC/Tasks/bqsr.wdl" as BQSR
import "MayomicsVC/src/wdl/gatk/HaplotyperVC/Tasks/haplotyper.wdl" as HAPLOTYPER

workflow CallHaplotyperVCTasks {

   Array[String] GenomicIntervals

   scatter (interval in GenomicIntervals) {
      call BQSR.bqsrTask as bqsr {
         input:
            GenomicInterval = interval
      }
      
      call HAPLOTYPER.variantCallingTask as haplotype { 
         input:
            InputBams = bqsr.OutputBams,
            InputBais = bqsr.OutputBais,
            GenomicInterval = interval
      }
  }
      
}
