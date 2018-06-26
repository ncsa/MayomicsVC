###########################################################################################


task ReadMappingTask {
        
   File Input_Read1		# Input Read File		 (REQUIRED)
   File Input_Read2		# Input Read File		 (Optional)
   String sampleName		# Name of the Sample
   File RefFasta		# Reference FASTA file		 (REQUIRED)

   File Ref_Amb_File		#
   File Ref_Dict_File		#
   File Ref_Ann_File		#
   File Ref_Bwt_File		# These are reference files that are provided as implicit inputs
   File Ref_Fai_File		# to the WDL Tool to help perform the alignment
   File Ref_Pac_File		#
   File Ref_Sa_File		#

   String BWA			# Variable path to BWA MEM Tool
   String SAMTOOL               # variable path to Samtools
   
   command <<<

      ${BWA} mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} | ${SAMTOOL} view -@ 17 -bSu -> ${sampleName}.aligned.bam
  
    >>>
   
   # The output block is where the output of the program is stored.
   output {
      File Aligned_Bam = "${sampleName}.aligned.bam"
      String Global_sampleName = "${sampleName}"
      Array[String] Collect_AlignedBam = [Global_sampleName, Aligned_Bam]
   }

   runtime {
      
      continueOnReturnCode: true
   }

}



workflow CallReadMappingTask {

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

   File InputSamplesFile

   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   scatter(sample in inputsamples) {
   
      call ReadMappingTask {
         input :
            sampleName = sample[0],
            Input_Read1 = sample[1],
            Input_Read2 = sample[2],

            RefFasta = RefFasta,
            Ref_Amb_File = Ref_Amb_File,
            Ref_Dict_File = Ref_Dict_File,
            Ref_Ann_File = Ref_Ann_File,
            Ref_Bwt_File = Ref_Bwt_File,
            Ref_Fai_File = Ref_Fai_File,
            Ref_Pac_File = Ref_Pac_File,
            Ref_Sa_File = Ref_Sa_File,

            BWA = BWA,
            SAMTOOL = SAMTOOL 
      }
   
   }

   output {
      Array[Array[String]] Global_AlignedBam = ReadMappingTask.Collect_AlignedBam
   }
   
}
