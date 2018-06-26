###########################################################################################

##              This WDL script performs alignment using BWA Mem                ##

##                              Script Options
#       -t        "Number of Threads"                         (Optional)
#       -P        "Single Ended Reads specification"          (Required)
#       -l        "Left Fastq File"                           (Required)
#       -r        "Right Fastq File"                          (Optional)
#       -G        "Reference Genome"                          (Required)
#       -s        "Name of the sample"                        (Optional)
#       -S        "Path to the Sentieon Tool"                 (Required)
#       -L        "Sentieon License File"                     (Required)
#       -g        "Group"                                     (Required)
#       -p        "Platform"                                  (Required)

###########################################################################################

task alignmentTask {

   File RefFasta                   # Reference Input Fasta File
   File InputRead1                 # Input Read File           
   String InputRead2               # Input Read File           
   String SampleName               # Name of the Sample

   File RefAmb                     
   File RefAnn                     
   File RefBwt                     # These are reference files that are provided as implicit inputs
   File RefPac                     # to the WDL Tool to help perform the alignment
   File RefSa              

   String SentieonLicense          # Sentieon License server
   String Sentieon                 # Path to Sentieon

   String Group
   String Platform
   String DebugMode                # Flag to enable Debug Mode
   String Threads                  # Specifies the number of thread required per run

   Boolean PairedEnd               # Variable to check if single ended or not

   File AlignmentScript            # Bash script which is called inside the WDL script

   command {

      /bin/bash ${AlignmentScript} -L ${SentieonLicense} -P ${PairedEnd} -g ${Group} -l ${InputRead1} -r ${InputRead2} -s ${SampleName} -p ${Platform} -G ${RefFasta} -S ${Sentieon} -t ${Threads} ${DebugMode}

   }

   output {

      File AlignedSortedBam = "${SampleName}.aligned.sorted.bam"
      File AlignedSortedBamIdx = "${SampleName}.aligned.sorted.bam.bai"

   }

} 

