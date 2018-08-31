#################################################################################################

##              This WDL script marks the duplicates on input sorted BAMs                ##

#                              Script Options
#      -b        "Input BAM File"                            (Required)      
#      -s        "Name of the sample"                        (Optional)
#      -t        "Number of Threads"                         (Optional)
#      -S        "Path to the Sentieon Tool"                 (Required)
#      -O        "Directory for the Output"                  (Required)
#      -L        "Sentieon License File"                     (Required)
#      -e        "Path to the environmental profile          (Required)
#      -d        "Debug Mode Toggle"                         (Optional)

#################################################################################################

task dedupTask {

   File InputAlignedSortedBam             # Input Sorted BAM File
   File InputAlignedSortedBamBai          # Input Sorted Bam Index File

   String SampleName                      # Name of the Sample

   String Sentieon                        # Variable path to Sentieon 
   String SentieonLicense                 # License Server Information

   String SentieonThreads                 # Specifies the number of thread required per run
   String DebugMode                       # Variable to check whether Debud Mode is on

   File DedupScript                       # Bash script that is called inside the WDL script
   File EnvProfile                        # File containing the environmental profile variables

   command {

      /bin/bash ${DedupScript} -L ${SentieonLicense} -b ${InputAlignedSortedBam} -s ${SampleName} -S ${Sentieon} -t ${SentieonThreads} -e ${EnvProfile} ${DebugMode}

   }

   output {

      File AlignedSortedDeduppedBam = "${SampleName}.aligned.sorted.deduped.bam"
      File AlignedSortedDeduppedBamBai = "${SampleName}.aligned.sorted.deduped.bam.bai"

   }
}
