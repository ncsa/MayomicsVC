#########################################################################################################

###       This WDL script performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "../Tasks/alignment.wdl" as ALIGN

workflow CallalignmentTask {

   Array[Array[File]] InputReads   # One lane per subarray with one or two input reads
   String SampleName               # Name of the Sample
   String Group                    # starting read group string
   String Platform                 # sequencing platform for read group
   Boolean PairedEnd               # Variable to check if single ended or not

   File Ref                        # Reference Genome
   File RefAmb                     # reference file index
   File RefAnn                     # reference file index
   File RefBwt                     # reference file index
   File RefPac                     # reference file index
   File RefSa                      # reference file index

   String Sentieon                 # Path to Sentieon
   String SentieonThreads          # Specifies the number of thread required per run

   File AlignmentScript            # Bash script which is called inside the WDL script
   File AlignEnvProfile            # File containing the environmental profile variables
   String ChunkSizeInBases         # The -K option for BWA MEM


   String DebugMode                # Flag to enable Debug Mode
   

   scatter (lane in InputReads) {

      if(PairedEnd) {
         call ALIGN.alignmentTask as ALIGN_paired {
            input:
               InputRead1=lane[0],
               InputRead2=lane[1],
               Ref=Ref,
               RefAmb=RefAmb,
               RefAnn=RefAnn,
               RefBwt=RefBwt,
               RefPac=RefPac,
               RefSa=RefSa,
               SampleName=SampleName,
               Group=Group,
               Platform=Platform,
               PairedEnd=PairedEnd,
               Sentieon=Sentieon,
               SentieonThreads=SentieonThreads,
               AlignmentScript=AlignmentScript,
               AlignEnvProfile=AlignEnvProfile,
               ChunkSizeInBases=ChunkSizeInBases,
               DebugMode=DebugMode
         }
      }

      if(!PairedEnd) {
         call ALIGN.alignmentTask as ALIGN_single {
            input:
               InputRead1=lane[0],
               InputRead2="null",
               Ref=Ref,
               RefAmb=RefAmb,
               RefAnn=RefAnn,
               RefBwt=RefBwt,
               RefPac=RefPac,
               RefSa=RefSa,
               SampleName=SampleName,
               Group=Group,
               Platform=Platform,
               PairedEnd=PairedEnd,
               Sentieon=Sentieon,
               SentieonThreads=SentieonThreads,
               AlignmentScript=AlignmentScript,
               AlignEnvProfile=AlignEnvProfile,
               ChunkSizeInBases=ChunkSizeInBases,
               DebugMode=DebugMode
         }
      }
   }

   output {
      # Unify outputs from scatter and filter out null entries 
      Array[File] AlignedSortedBams = select_all(flatten([ALIGN_paired.AlignedSortedBam,ALIGN_single.AlignedSortedBam]))
      Array[File] AlignedSortedBamBais = select_all(flatten([ALIGN_paired.AlignedSortedBamBai,ALIGN_single.AlignedSortedBamBai]))
   }

}
