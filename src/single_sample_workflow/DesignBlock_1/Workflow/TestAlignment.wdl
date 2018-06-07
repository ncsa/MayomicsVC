##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "DesignBlock_1/Tasks/CutAdaptTrimming.wdl" as CUTADAPTTRIM
import "DesignBlock_1/Tasks/Sentieon/Alignment.wdl" as ALIGNMENT
import "DesignBlock_1/Tasks/Sentieon/MarkDuplicates.wdl" as DEDUP 
import "DesignBlock_1/Tasks/Sentieon/Realignment.wdl" as REALIGN

workflow CallAlignmentStageTasks {
   
   call CUTADAPTTRIM.TrimInputSequencesTask 
    
   call ALIGNMENT.ReadMappingTask {
      input:
         sampleName = TrimInputSequencesTask.sampleName,
         Debug_Mode_EN = TrimInputSequencesTask.Debug_Mode_EN,
         Error_Logs = TrimInputSequencesTask.Error_Logs,
         Threads = TrimInputSequencesTask.Threads,
         Is_Single_End = TrimInputSequencesTask.Is_Single_End,
         OutDir = TrimInputSequencesTask.OutDir,
         Input_Read1 = TrimInputSequencesTask.TrimmedInputRead1,
         Input_Read2 = TrimInputSequencesTask.TrimmedInputRead2
   }
   
   call DEDUP.MarkDuplicatesTask {
      input:
         Sentieon = ReadMappingTask.Sentieon,
         sampleName = ReadMappingTask.sampleName,
         Debug_Mode_EN = ReadMappingTask.Debug_Mode_EN,
         Error_Logs = ReadMappingTask.Error_Logs,
         OutDir = ReadMappingTask.OutDir,
         Threads = ReadMappingTask.Threads,
         InputBam  = ReadMappingTask.SortBam
   }

   call REALIGN.RealignmentTask {
      input:
         Sentieon = MarkDuplicatesTask.Sentieon,
         sampleName = MarkDuplicatesTask.sampleName,
         Debug_Mode_EN = MarkDuplicatesTask.Debug_Mode_EN,
         InputBam = MarkDuplicatesTask.AlignedSortedDeduppedBam,
         Error_Logs = MarkDuplicatesTask.Error_Logs,
         OutDir = MarkDuplicatesTask.OutDir,
         Threads = MarkDuplicatesTask.Threads,
         Ref_Amb_File = ReadMappingTask.Ref_Amb_File,
         Ref_Dict_File = ReadMappingTask.Ref_Dict_File,
         Ref_Ann_File = ReadMappingTask.Ref_Ann_File,
         Ref_Bwt_File = ReadMappingTask.Ref_Bwt_File,
         Ref_Fai_File = ReadMappingTask.Ref_Fai_File,
         Ref_Pac_File = ReadMappingTask.Ref_Pac_File,
         Ref_Sa_File = ReadMappingTask.Ref_Sa_File
   }   
   
}
