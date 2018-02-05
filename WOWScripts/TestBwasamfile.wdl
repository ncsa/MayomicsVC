# The Task block is where the variables and the functions are defined for performing a certain task

task ReadMappingTask {

   File Input_Read1             # Input Read File                (REQUIRED)
   File Input_Read2             # Input Read File                (Optional)
   String sampleName            # Name of the Sample
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
   String dollar = "$"          # Variable to access internal bash variables

   command <<<

      # BWA Mem Tool is used to create aligned SAM file from the input FASTA File
      # PIPESTATUS is an internal bash variable which holds the exit code of commands in the pipe
      ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} | ${SAMTOOL} view -@ 17 -bSu -> ${sampleName}.aligned.bam; B=(${dollar}{PIPESTATUS[*]})

   >>>

   # The output block is where the output of the program is stored.
   output {
      File Aligned_Bam = "${sampleName}.aligned.bam"
   }

   # Runtime block specifies the Cromwell Engine of runtime attributes to customize the environment for the call
   runtime {
      # Even if the command in task has a non zero exit code continue with the other tasks
      continueOnReturnCode: true
      failOnStderr: false
   }

}  # End of task block

workflow CallReadMappingTask{

   File Input_Read1            
   File Input_Read2            
   String sampleName           
   File RefFasta              

   File Ref_Amb_File            
   File Ref_Dict_File           
   File Ref_Ann_File            
   File Ref_Bwt_File           
   File Ref_Fai_File           
   File Ref_Pac_File           
   File Ref_Sa_File          

   String BWA
   String SAMTOOL
  
   call ReadMappingTask {
      input :

         RefFasta = RefFasta,
         Ref_Amb_File = Ref_Amb_File,
         Ref_Dict_File = Ref_Dict_File,
         Ref_Ann_File = Ref_Ann_File,
         Ref_Bwt_File = Ref_Bwt_File,
         Ref_Fai_File = Ref_Fai_File,
         Ref_Pac_File = Ref_Pac_File,
         Ref_Sa_File = Ref_Sa_File,

         Input_Read1 = Input_Read1,
         Input_Read2 = Input_Read2,
         RefFasta = RefFasta,
         sampleName = sampleName,

         BWA = BWA,
         SAMTOOL = SAMTOOL

   }

   output {
      File Global_Aligned_Bam = ReadMappingTask.Aligned_Bam
   }

}
