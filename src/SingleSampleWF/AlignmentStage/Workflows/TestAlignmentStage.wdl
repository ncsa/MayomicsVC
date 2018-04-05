#########################################################################################################
###              This WDL script is used to run the Alignment steps as individual modules              ##
#########################################################################################################

import "AlignmentStage_WDL/Tasks/PreExecQC.wdl" as PREEXECQC
import "AlignmentStage_WDL/Tasks/BWAMemSamtoolView.wdl" as BWAMEMSAMTOOLVIEW
import "AlignmentStage_WDL/Tasks/Novosort.wdl" as NOVOSORT
import "AlignmentStage_WDL/Tasks/PicardMarkDuplicates.wdl" as PICARDMARKDUPLICATES
import "AlignmentStage_WDL/Tasks/EndOfBlockNotify.wdl" as ENDOFBLOCKNOTIFY


workflow CallAlignmentStageTasks {

   call PREEXECQC.QualityControlTask {
      input :
         failure_logs = Failure_Logs
      }

 
   call BWAMEMSAMTOOLVIEW.ReadMappingTask {
      input :
         RefFasta = QualityControlTask.RefFasta,
         Ref_Amb_File = QualityControlTask.Ref_Amb_File,
         Ref_Dict_File = QualityControlTask.Ref_Dict_File, 
         Ref_Ann_File = QualityControlTask.Ref_Ann_File,            
         Ref_Bwt_File = QualityControlTask.Ref_Bwt_File,            
         Ref_Fai_File = QualityControlTask.Ref_Fai_File,            
         Ref_Pac_File = QualityControlTask.Ref_Pac_File,            
         Ref_Sa_File = QualityControlTask.Ref_Sa_File,           
         sampleName = sample[0],         
         Input_Read1 = sample[1],         
         Input_Read2 = sample[2],
         BWA = QualityControlTask.BWA,
         SAMTOOL = QualityControlTask.SAMTOOL,
         Exit_Code = Capture_Exit_Code,
         Failure_Logs = QualityControlTask.Failure_Logs,
         
         # a dummy variable to enforce sequentiality b/w pre-exec QC and alignment 
         # due to the absence of explicit data dependency
         DummyVar = QualityControlTask.DummyVar
   }   
    
   call NOVOSORT.NovosortTask {
      input :
         sampleName = sample[0],
         SORT = QualityControlTask.SORT,
         Aligned_Bam = ReadMappingTask.Aligned_Bam,
         Exit_Code = Capture_Exit_Code,
         Failure_Logs = Failure_Logs
   }

   call PICARDMARKDUPLICATES.MarkDuplicatesTask {
      input :
         sampleName = sample[0],
         JAVA = QualityControlTask.JAVA,
         PICARD = QualityControlTask.PICARD,
         Aligned_Sorted_Bam = NovosortTask.Aligned_Sorted_Bam,
         Exit_Code = Capture_Exit_Code,
         Failure_Logs = Failure_Logs
   } 



   call ENDOFBLOCKNOTIFY.EndOfStageEmailTask {
      input :
         Email = MarkDuplicatesTask.Notify_EndofAlignment,
         Failure_Logs = Failure_Logs
   }

} # End of Workflow block


