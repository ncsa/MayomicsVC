############

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentStage_WDL/TestBwsamfile.wdl" as BWA

task NovosortTask {
   Array[File] Aligned_Bam                        # Input BAM File
   String sampleName                       # Name of the Sample
   String Exit_Code                        # File to capture exit code
   String SORT                             # Variable path to Novosort
   String Failure_Logs                     # Variable to capture Failure Reports


   command {

      # Check to see if input files are non-zero
      [ -s ${Aligned_Bam} ] || echo "Input BAM File is Empty" >> ${Failure_Logs}

      # Novosort Tools is used to created sort BAM Files 
      ${SORT} -c 36 -i -o ${sampleName}.aligned.sorted.bam ${sep=',' Aligned_Bam}

      # The 'if' check to see if any of the samples have failed this step          
      if [ $? -ne 0 ]; then
         echo '${sampleName} has failed at the Novosort Step' >> ${Exit_Code}
      fi

      [ ! -f ${sampleName}.aligned.sorted.bam ] && echo "aligned sorted bam not created" >> ${Failure_Logs}
   }

   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs   
   output {
      File Aligned_Sorted_Bam = "${sampleName}.aligned.sorted.bam"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: true
   }
} # End of Task Block


workflow CallNovosortTask {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   #The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)


################### BIOLERPLATE FOR BWASAM ######################################################################

   File Input_Read1             # Input Read File                (REQUIRED)
   File Input_Read2             # Input Read File                (Optional)
   File RefFasta                # Reference FASTA file           (REQUIRED)

   File Ref_Amb_File            #
   File Ref_Dict_File           #
   File Ref_Ann_File            #
   File Ref_Bwt_File            # These are reference files that are provided as implicit inputs
   File Ref_Fai_File            # to the WDL Tool to help perform the alignment
   File Ref_Pac_File            #
   File Ref_Sa_File             #

   String BWA                   # Variable path to BWA MEM Tool
   String SAMTOOL               # variable path to Samtools

################################################################################################################


   call BWA.CallReadMappingTask {

      input :   
         RefFasta = ReadMappingTask.RefFasta,
         Ref_Amb_File = ReadMappingTask.Ref_Amb_File,
         Ref_Dict_File = ReadMappingTask.Ref_Dict_File,
         Ref_Ann_File = ReadMappingTask.Ref_Ann_File,
         Ref_Bwt_File = ReadMappingTask.Ref_Bwt_File,
         Ref_Fai_File = ReadMappingTask.Ref_Fai_File,
         Ref_Pac_File = ReadMappingTask.Ref_Pac_File,
         Ref_Sa_File = ReadMappingTask.Ref_Sa_File

   }

   # The scatter function performs operations on multiple samples in parallel
   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {

      # Novosort is included as a sub task and it is called inside the workflow
      call NovosortTask {
         input :
            sampleName = sample[0],
            Aligned_Bam = CallReadMappingTask.Global_Aligned_Bam
      }
  }
}
