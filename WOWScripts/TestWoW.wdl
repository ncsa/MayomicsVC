##################################################################################################################

import "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/AlignmentStage_WDL/TestBwasamfile.wdl" as BWA

task NovosortTask {
   File Aligned_Bam                 # Input BAM File
   String sampleName                       # Name of the Sample
   String Exit_Code                        # File to capture exit code
   String SORT                             # Variable path to Novosort
   String Failure_Logs                     # Variable to capture Failure Reports


   command {


      # Novosort Tools is used to created sort BAM Files 
      ${SORT} -c 18 -i -o ${sampleName}.aligned.sorted.bam ${Aligned_Bam}

   }

   # The output block is where the output of the program is stored.
   # glob function is used to capture the multi sample outputs   
   output {
      File Aligned_Sorted_Bam = "${sampleName}.aligned.sorted.bam"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      continueOnReturnCode: false
   }
} # End of Task Block


workflow CallNovosortTask {
   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   #The 2-D array stores information of all the samples 
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   # The scatter operation has an implied gather step at the end
   scatter(sample in inputsamples) {

      call BWA.CallReadMappingTask as sub_CallReadMappingTask {
         input :
            sampleName = sample[0],
            Input_Read1 = sample[1],
            Input_Read2 = sample[2]
      }

      # Novosort is included as a sub task and it is called inside the workflow
      call NovosortTask {
         input :
            sampleName = sample[0],
            Aligned_Bam = sub_CallReadMappingTask.Global_Aligned_Bam
      }
   }
   
}
