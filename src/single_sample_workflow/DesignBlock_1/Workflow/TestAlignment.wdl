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
         sampleName = TrimInputSequencesTask.sName,
         Debug_Mode_EN = TrimInputSequencesTask.Debug_Mode,
         Error_Logs = TrimInputSequencesTask.ErrLogs,
         Threads = TrimInputSequencesTask.ThreadCount,
         Is_Single_End = TrimInputSequencesTask.SingleEnd,
         OutDir = TrimInputSequencesTask.OutputDir,
         Input_Read1 = TrimInputSequencesTask.TrimmedInputRead1,
         Input_Read2 = TrimInputSequencesTask.TrimmedInputRead2
   }
   
   call DEDUP.MarkDuplicatesTask {
      input:
         Sentieon_License = ReadMappingTask.LicenseFile,
         Sentieon = ReadMappingTask.SentieonPath,
         sampleName = ReadMappingTask.sName,
         Debug_Mode_EN = ReadMappingTask.DebugMode,
         Error_Logs = ReadMappingTask.ErrLogs,
         OutDir = ReadMappingTask.OutputDir,
         Threads = ReadMappingTask.ThreadCount,
         InputBam  = ReadMappingTask.SortBam,
         InputBamIdx = ReadMappingTask.SortBamIdx
   }

   call REALIGN.RealignmentTask {
      input:
         Sentieon_License = MarkDuplicatesTask.LicenseFile,
         Sentieon = MarkDuplicatesTask.SentieonPath,
         sampleName = MarkDuplicatesTask.sName,
         Debug_Mode_EN = MarkDuplicatesTask.DebugMode,
         InputBam = MarkDuplicatesTask.AlignedSortedDeduppedBam,
         InputBamIdx = MarkDuplicatesTask.AlignedSortedDeduppedBamIndex,
         Error_Logs = MarkDuplicatesTask.ErrLogs,
         OutDir = MarkDuplicatesTask.OutputDir,
         Threads = MarkDuplicatesTask.ThreadCount,
         RefFasta = ReadMappingTask.FastaRef,
         Ref_Amb_File = ReadMappingTask.RefAmbFile,
         Ref_Dict_File = ReadMappingTask.RefDictFile,
         Ref_Ann_File = ReadMappingTask.RefAnnFile,
         Ref_Bwt_File = ReadMappingTask.RefBwtFile,
         Ref_Fai_File = ReadMappingTask.RefFaiFile,
         Ref_Pac_File = ReadMappingTask.RefPacFile,
         Ref_Sa_File = ReadMappingTask.RefSaFile
   }   
   
}
