########################################################################################################

## This wdl script checks if the inputs and executables are present and are non-zero ##

########################################################################################################

task QualityControlTask {

   # Variables that hold the path to executables and inputs

   #Executables 
   String bwa
   String samtool
   String Java
   String sort
   String picard
   String gatk
   Int Flag = 0

   #Inputs
   File ref_fasta

   #Reference Files 
   File Ref_amb_File
   File Ref_dict_File          
   File Ref_ann_File            
   File Ref_bwt_File            
   File Ref_fai_File            
   File Ref_pac_File            
   File Ref_sa_File
   #Variable to capture Failure reports
   String failure_logs


   command <<<
   
      # To check if the executables and inputs are present
      [ ! -f ${bwa} ] && echo "BWA does not exist" >> ${failure_logs}
      [ ! -f ${samtool} ] && echo "Samtools does not exist" >> ${failure_logs}
      [ ! -f ${Java} ] && echo "Java does not exist" >> ${failure_logs}
      [ ! -f ${sort} ] && echo "Novosort does not exist" >> ${failure_logs}
      [ ! -f ${picard} ] && echo "Picard does not exist" >> ${failure_logs}
      [ ! -f ${gatk} ] && echo "GATK does not exist" >> ${failure_logs}

      [ ! -f ${ref_fasta} ] && echo "Reference Fasta File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_amb_File} ] && echo "Reference Amb File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_dict_File} ] && echo "Reference Dict File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_ann_File} ] && echo "Reference Ann File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_bwt_File} ] && echo "Reference Bwt File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_fai_File} ] && echo "Reference Fai File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_pac_File} ] && echo "Reference Pac File does does not exist" >> ${failure_logs}
      [ ! -f ${Ref_sa_File} ] && echo "Reference Sa File does does not exist" >> ${failure_logs}  

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
      String GATK = gatk
      String Failure_Logs = failure_logs
      File RefFasta = ref_fasta
      File Ref_Amb_File = Ref_amb_File
      File Ref_Dict_File = Ref_dict_File 
      File Ref_Ann_File = Ref_ann_File
      File Ref_Bwt_File = Ref_bwt_File
      File Ref_Fai_File = Ref_fai_File
      File Ref_Pac_File = Ref_pac_File
      File Ref_Sa_File = Ref_sa_File
   }

}
