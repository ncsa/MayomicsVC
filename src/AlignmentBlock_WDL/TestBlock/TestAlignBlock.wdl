#########################################################################################################
###              This WDL script is used to run the Alignment steps as individual modules              ##
#########################################################################################################

### ??? figure out how to convert from variable to an explicit path in the import command ??? ###
import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/PreExec_QC.wdl" as PreQC
import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/Bwa_Sam.wdl" as BWA
import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/Novosort.wdl" as NSORT
import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/PicardMD.wdl" as PICARD
import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentBlock_WDL/Tasks/EndofBlock_Notify.wdl" as EMAIL



workflow AlignBlock_Run {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   # The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   String Capture_Exit_Code
   String Failure_Logs
   
   # Task to check if the executables and inputs exists and are non-zero
   call PreQC.PreExec_QC {
      input :
         failure_logs = Failure_Logs
      }

   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {
 
      # BWA Mem is included as a sub task and it is called inside the workflow
      call BWA.BWA_Mem {
         input :
            RefFasta = PreExec_QC.RefFasta,
            Ref_Amb_File = PreExec_QC.Ref_Amb_File,
            Ref_Dict_File = PreExec_QC.Ref_Dict_File, 
            Ref_Ann_File = PreExec_QC.Ref_Ann_File,            
            Ref_Bwt_File = PreExec_QC.Ref_Bwt_File,            
            Ref_Fai_File = PreExec_QC.Ref_Fai_File,            
            Ref_Pac_File = PreExec_QC.Ref_Pac_File,            
            Ref_Sa_File = PreExec_QC.Ref_Sa_File,           
            sampleName = sample[0],         
            Input_Read1 = sample[1],         
            Input_Read2 = sample[2],
            BWA = PreExec_QC.BWA,
            SAMTOOL = PreExec_QC.SAMTOOL,
            Exit_Code = Capture_Exit_Code,
            Failure_Logs = PreExec_QC.Failure_Logs,
            
            # a dummy variable to enforce sequentiality b/w pre-exec QC and alignment 
            # due to the absence of explicit data dependency
            DummyVar = PreExec_QC.DummyVar
      }   
    
      call NSORT.Novosort {
         input :
            sampleName = sample[0],
            SORT = PreExec_QC.SORT,
            Aligned_Bam = BWA_Mem.Aligned_Bam,
            Exit_Code = Capture_Exit_Code,
            Failure_Logs = Failure_Logs
      }

      call PICARD.Picard_MarkDuplicates {
         input :
            sampleName = sample[0],
            JAVA = PreExec_QC.JAVA,
            PICARD = PreExec_QC.PICARD,
            Aligned_Sorted_Bam = Novosort.Aligned_Sorted_Bam,
            Exit_Code = Capture_Exit_Code,
            Failure_Logs = Failure_Logs
      }      

   } # End of scatter block

   call EMAIL.EndofBlock_Notify {
      input :
         Email = Picard_MarkDuplicates.Notify_EndofAlignment,
         Failure_Logs = Failure_Logs
   }

} # End of Workflow block


