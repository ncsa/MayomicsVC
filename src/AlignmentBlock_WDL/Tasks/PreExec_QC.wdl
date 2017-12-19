########################################################################################################

## This wdl script checks if the inputs and executables are present and are non-zero ##

########################################################################################################

task PreExec_QC {

   # Variables that hold the path to executables and inputs

   #Executables 
   String bwa
   String samtool
   String Java
   String sort
   String picard
   Int Flag = 0

   #Inputs
   File ref_fasta

   #Variable to capture Failure reports
   String failure_logs


   command <<<
   
      # To check if the executables and inputs are present
      [ ! -f ${bwa} ] && echo "BWA does not exist" >> ${failure_logs}
      [ ! -f ${samtool} ] && echo "Samtools does not exist" >> ${failure_logs}
      [ ! -f ${Java} ] && echo "Java does not exist" >> ${failure_logs}
      [ ! -f ${sort} ] && echo "Novosort does not exist" >> ${failure_logs}
      [ ! -f ${picard} ] && echo "Picard does not exist" >> ${failure_logs}
      [ ! -f ${ref_fasta} ] && echo "Reference File does does not exist" >> ${failure_logs}
   
      # To check if the Input Files are non-zero
      [ -s ${ref_fasta} ] || echo "Reference File is Empty" >> ${failure_logs}
   
   >>>


   output {

      Int DummyVar = Flag
      String BWA = bwa
      String SAMTOOL = samtool
      String JAVA = Java
      String SORT = sort
      String PICARD = picard
      String Failure_Logs = failure_logs
      File RefFasta = ref_fasta
   }

}
