########

task ReadMappingTask {

   String sampleName

   File RefFasta                   # Reference FASTA file        (REQUIRED)
   File Ref_Amb_File               #
   File Ref_Ann_File               #
   File Ref_Bwt_File               # These are reference files that are provided as implicit inputs
   File Ref_Fai_File               # to the WDL Tool to help perform the alignment
   File Ref_Pac_File               #
   File Ref_Sa_File                #

   File Input_Read1                # Input Read File             (REQUIRED)
   File Input_Read2                # Input Read File             (Optional)

   String BWA                      # Variable path to BWA MEM Tool
   String SAMTOOLS                 # variable path to Samtools
   String BashScriptPath
   String Failure_Logs             # Variable to capture Failure reports

   command <<<
      /bin/bash ${BashScriptPath} '${Failure_Logs}' '${BWA}' '${RefFasta}' '${Input_Read1}' '${Input_Read2}' '${SAMTOOLS}' '${sampleName}'
   >>>

   output {
      File Aligned_Bam = "${sampleName}.aligned.bam"
   }
}






workflow CallBWAMemSamtoolsViewTask {

   # The InputSamplesFile is a variable that stores information on various samples
   File InputSamplesFile

   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

   scatter(sample in inputsamples) {

      call ReadMappingTask {
         input :
            sampleName = sample[0],
            Input_Read1 = sample[1],
            Input_Read2 = sample[2],
         
      }
   }
}
