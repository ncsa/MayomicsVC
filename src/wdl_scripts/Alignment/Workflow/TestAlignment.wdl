##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl" as CUTADAPTTRIM
import "src/wdl_scripts/Alignment/TestTasks/Runalignment.wdl" as ALIGNMENT
import "src/wdl_scripts/Alignment/Tasks/dedup.wdl" as DEDUP 

workflow CallAlignmentTasks {
   
   Boolean Trimming
   Array[Array[File]] InputReads

   if(Trimming) {

      call CUTADAPTTRIM.RunTrimSequencesTask as trimseq {
         input:
            InputReads = InputReads
      } 
       
      call ALIGNMENT.RunAlignmentTask as align_w_trim {
         input:
            InputReads = trimseq.Outputs
      }
   }
   if(!Trimming) {
      
      call ALIGNMENT.RunAlignmentTask as align_wo_trim {
         input:
            InputReads = InputReads
      }
   }
   
   call DEDUP.dedupTask as dedup {
      input:
         InputBams = select_all([align_w_trim.OutputBams,align_wo_trim.OutputBams]),
         InputBais = select_all([align_w_trim.OutputBais,align_wo_trim.OutputBais])
   }

   output {
     
      File InputBams = dedup.OutputBams
      File InputBais = dedup.OutputBais
   }    
    
}
