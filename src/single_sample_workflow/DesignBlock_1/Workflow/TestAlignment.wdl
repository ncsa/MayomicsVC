##########################################################################################################
####              This WDL script is used to run the Alignment steps as individual modules              ##
##########################################################################################################

import "DesignBlock_1/Tasks/trim_sequences.wdl" as CUTADAPTTRIM
import "DesignBlock_1/Tasks/Sentieon/alignment.wdl" as ALIGNMENT
import "DesignBlock_1/Tasks/Sentieon/dedup.wdl" as DEDUP 
import "DesignBlock_1/Tasks/Sentieon/realignment.wdl" as REALIGN

workflow CallAlignmentStageTasks {
   
   call CUTADAPTTRIM.trimsequencesTask 
    
   call ALIGNMENT.alignmentTask {
      input:
         sampleName = trimsequencesTask.sName,
         Debug_Mode_EN = trimsequencesTask.Debug_Mode,
         Error_Logs = trimsequencesTask.ErrLogs,
         Threads = trimsequencesTask.ThreadCount,
         Is_Single_End = trimsequencesTask.SingleEnd,
         OutDir = trimsequencesTask.OutputDir,
         Input_Read1 = trimsequencesTask.TrimmedInputRead1,
         Input_Read2 = trimsequencesTask.TrimmedInputRead2
   }
   
   call DEDUP.dedupTask {
      input:
         Sentieon_License = alignmentTask.LicenseFile,
         Sentieon = alignmentTask.SentieonPath,
         sampleName = alignmentTask.sName,
         Debug_Mode_EN = alignmentTask.DebugMode,
         Error_Logs = alignmentTask.ErrLogs,
         OutDir = alignmentTask.OutputDir,
         Threads = alignmentTask.ThreadCount,
         InputBam  = alignmentTask.SortBam,
         InputBamIdx = alignmentTask.SortBamIdx
   }

   call REALIGN.realignmentTask {
      input:
         Sentieon_License = dedupTask.LicenseFile,
         Sentieon = dedupTask.SentieonPath,
         sampleName = dedupTask.sName,
         Debug_Mode_EN = dedupTask.DebugMode,
         InputBam = dedupTask.AlignedSortedDeduppedBam,
         InputBamIdx = dedupTask.AlignedSortedDeduppedBamIndex,
         Error_Logs = dedupTask.ErrLogs,
         OutDir = dedupTask.OutputDir,
         Threads = dedupTask.ThreadCount,
         RefFasta = alignmentTask.FastaRef,
         Ref_Amb_File = alignmentTask.RefAmbFile,
         Ref_Dict_File = alignmentTask.RefDictFile,
         Ref_Ann_File = alignmentTask.RefAnnFile,
         Ref_Bwt_File = alignmentTask.RefBwtFile,
         Ref_Fai_File = alignmentTask.RefFaiFile,
         Ref_Pac_File = alignmentTask.RefPacFile,
         Ref_Sa_File = alignmentTask.RefSaFile
   }   
   
}
