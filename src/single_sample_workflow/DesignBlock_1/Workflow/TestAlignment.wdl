##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "DesignBlock_1/Tasks/CutAdaptTrimming.wdl" as CUTADAPTTRIM
import "DesignBlock_1/Tasks/Sentieon/BWAMemSamtoolsView_Sentieon.wdl" as BWAMEMSAMTOOLSVIEW
import "DesignBlock_1/Tasks/Sentieon/MarkDuplicates_Sentieon.wdl" as DEDUP 


workflow CallAlignmentStageTasks {

   call CUTADAPTTRIM.TrimInputSequencesTask 

   if (!TrimInputSequencesTask.Is_Single_End) {
         call BWAMEMSAMTOOLSVIEW.ReadMappingTask {
         input:
            Input_Read1 = TrimInputSequencesTask.TrimmedInputRead1,
            Input_Read2 = TrimInputSequencesTask.TrimmedInputRead2
      }

      call DEDUP.MarkDuplicatesTask {
      input:
         InputBam  = ReadMappingTask.SortBam
   }

   }

   if (TrimInputSequencesTask.Is_Single_End) {
      call BWAMEMSAMTOOLSVIEW.ReadMappingTask {
         input:
            Input_Read1 = TrimInputSequencesTask.TrimmedInputRead
      }
   
   call DEDUP.MarkDuplicatesTask {
      input:
         InputBam  = ReadMappingTask.SortBam
    }

   }

}
