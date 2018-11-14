#########################################################################################################

###       This WDL script performs BWA to create sam files and converts to bam using Samtools       ##

#########################################################################################################

import "src/wdl_scripts/Alignment/Tasks/alignment.wdl" as ALIGN

workflow RunAlignmentTask {

   Array[Array[File]] InputReads   # One lane per subarray with one or two input reads
   String SampleName               # Name of the Sample
   String Platform                 # sequencing platform for read group
   String Library                  # sequencing library for read group
   String CenterName               # sequencing center name for read group
   Array[String] PlatformUnit      # One platform unit per alignment task
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
   String BWAExtraOptionsString    # String of extra options for BWA. This can be an empty string.


   String DebugMode                # Flag to enable Debug Mode
   
   Array[Int] Indexes = range(length(InputReads))

   scatter (idx in Indexes) {

      if(PairedEnd) {
         call ALIGN.alignmentTask as ALIGN_paired {
            input:
               InputRead1=InputReads[idx][0],
               InputRead2=InputReads[idx][1],
               Ref=Ref,
               RefAmb=RefAmb,
               RefAnn=RefAnn,
               RefBwt=RefBwt,
               RefPac=RefPac,
               RefSa=RefSa,
               SampleName=SampleName,
               Platform=Platform,
               Library=Library,
               CenterName=CenterName,
               PlatformUnit=PlatformUnit[idx],
               PairedEnd=PairedEnd,
               Sentieon=Sentieon,
               SentieonThreads=SentieonThreads,
               AlignmentScript=AlignmentScript,
               AlignEnvProfile=AlignEnvProfile,
               ChunkSizeInBases=ChunkSizeInBases,
               BWAExtraOptionsString=BWAExtraOptionsString,
               DebugMode=DebugMode
         }
      }

      if(!PairedEnd) {
         call ALIGN.alignmentTask as ALIGN_single {
            input:
               InputRead1=InputReads[idx][0],
               InputRead2="null",
               Ref=Ref,
               RefAmb=RefAmb,
               RefAnn=RefAnn,
               RefBwt=RefBwt,
               RefPac=RefPac,
               RefSa=RefSa,
               SampleName=SampleName,
               Platform=Platform,
               Library=Library,
               CenterName=CenterName,
               PlatformUnit=PlatformUnit[idx],
               PairedEnd=PairedEnd,
               Sentieon=Sentieon,
               SentieonThreads=SentieonThreads,
               AlignmentScript=AlignmentScript,
               AlignEnvProfile=AlignEnvProfile,
               ChunkSizeInBases=ChunkSizeInBases,
               BWAExtraOptionsString=BWAExtraOptionsString,
               DebugMode=DebugMode
         }
      }
   }

   output {
      # Unify outputs from scatter and filter out null entries 
      Array[File] OutputBams = select_all(flatten([ALIGN_paired.OutputBams,ALIGN_single.OutputBams]))
      Array[File] OutputBais = select_all(flatten([ALIGN_paired.OutputBais,ALIGN_single.OutputBais]))
   }

}
