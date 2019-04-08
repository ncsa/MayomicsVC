###############################################################################################
####              This WDL script is used to run Alignment and HaplotyperVC blocks together  ##
###############################################################################################

workflow GermlineMasterWF {

   Boolean Trimming
   Boolean MarkDuplicates
   Array[Array[File]] NormalInputReads # One lane per subarray with one or two input reads
   Boolean PairedEnd               # Variable to check if single ended or not
   String SampleName               # Name of the Sample
   Array[String] PlatformUnit      # One platform unit per alignment task

   if(Trimming) {
        scatter (lane in NormalInputReads) {
           if(PairedEnd) {
              call trimsequencesTask as TRIMSEQ_paired {
                 input:
                    SampleName=SampleName,
                    InputRead1=lane[0],
                    InputRead2=lane[1]
              }
           }
           if(!PairedEnd) {
              call trimsequencesTask as TRIMSEQ_single {
                 input:
                    SampleName=SampleName,
                    InputRead1=lane[0],
                    InputRead2="null"
              }
           }
        }
        Array[Array[File]] trimOutputs = select_all(flatten([TRIMSEQ_paired.Outputs,TRIMSEQ_single.Outputs]))
   }
   
   Array[Array[File]] AlignInputReads = select_first([trimOutputs,NormalInputReads])
   Array[Int] Indexes = range(length(AlignInputReads))

   scatter (idx in Indexes) {
      if(PairedEnd) {
         call alignmentTask as ALIGN_paired {
            input:
               SampleName=SampleName,
               InputRead1=AlignInputReads[idx][0],
               InputRead2=AlignInputReads[idx][1],
               PlatformUnit=PlatformUnit[idx]
         }
      }
      if(!PairedEnd) {
         call alignmentTask as ALIGN_single {
            input:
               SampleName=SampleName,
               InputRead1=AlignInputReads[idx][0],
               InputRead2="null",
               PlatformUnit=PlatformUnit[idx]
         }
      }
   }

   Array[File] alignOutputBams = select_all(flatten([ALIGN_paired.OutputBams,ALIGN_single.OutputBams]))
   Array[File] alignOutputBais = select_all(flatten([ALIGN_paired.OutputBais,ALIGN_single.OutputBais]))

   call mergebamTask as merge {
      input:
         InputBams = alignOutputBams,
         InputBais = alignOutputBais
   }

   if(MarkDuplicates) {
      call dedupTask as dedup {
         input:
            InputBams = merge.OutputBams,
            InputBais = merge.OutputBais
      }
   }

   File DeliverAlignOutputBams = select_first([dedup.OutputBams,merge.OutputBams])
   File DeliverAlignOutputBais = select_first([dedup.OutputBais,merge.OutputBais])

   call deliverAlignmentTask as DAB {
      input:
         InputBams = DeliverAlignOutputBams,
         InputBais = DeliverAlignOutputBais
   }

   Array[String] GenomicIntervals

   scatter (interval in GenomicIntervals) {
      call bqsrTask as bqsr {
         input:
            InputBams = DeliverAlignOutputBams,
            InputBais = DeliverAlignOutputBais,
            GenomicInterval = interval
      }
      call variantCallingTask as haplotype {
         input:
            InputBams = bqsr.OutputBams,
            InputBais = bqsr.OutputBais,
            GenomicInterval = interval
      }
   }

   call mergegvcfsTask as mergegvcfs {
      input:
         InputGvcfs = haplotype.OutputVcf,
         InputIdxs = haplotype.OutputVcfIdx
   }

   call deliverHaplotyperVCTask as DHVC {
      input:
         InputVcf = mergegvcfs.OutputVcf,
         InputVcfIdx = mergegvcfs.OutputVcfIdx
    }
}

###############################################################################################
####                    The following parts define the WDL tasks needed                      ##
###############################################################################################


task trimsequencesTask {

   File InputRead1                 # Input Read File
   String InputRead2               # Input Read File

   String SampleName               # Name of the Sample

   File Adapters                   # Adapter FastA File

   String CutAdapt                 # Path to CutAdapt Tool
   String CutAdaptThreads          # Number of threads for cutadapt to use

   Boolean PairedEnd               # Variable to check if single ended or not

   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File TrimSeqScript              # Bash script which is called inside the WDL script

   String TrimSoftMemLimit         # Soft memory limit - nice shutdown
   String TrimHardMemLimit         # Hard memory limit - kill immediately

   String DebugMode                # Variable to check if Debug Mode is on or not


   command <<<
     source ${BashPreamble}
     /bin/bash ${TrimSeqScript} -P ${PairedEnd} -l ${InputRead1} -r ${InputRead2} -s ${SampleName} -A ${Adapters} -C ${CutAdapt} -t ${CutAdaptThreads} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${CutAdaptThreads}"
      s_vmem: "${TrimSoftMemLimit}"
      h_vmem: "${TrimHardMemLimit}"
   }

   output {
      Array[File] Outputs = glob("*.fastq.gz")
   }
}


task alignmentTask {
   String SampleName               # Name of the Sample
   String Platform                 # sequencing platform for read group
   String Library                  # Sequencing library for read group
   String PlatformUnit             # Platform unit / flowcell ID for read group
   String CenterName               # Name of the sequencing center for read group
   Boolean PairedEnd               # Variable to check if single ended or not
   File InputRead1                 # Input Read File
   String InputRead2               # Input Read File
   File Ref                        # Reference Genome
   File RefAmb                     # reference file index
   File RefAnn                     # reference file index
   File RefBwt                     # reference file index
   File RefPac                     # reference file index
   File RefSa                      # reference file index
   File BWAExe                     # Path to BWA executable
   String ChunkSizeInBases         # The -K option for BWA MEM
   String BWAExtraOptionsString    # String of extra options for BWA. This can be an empty string.
   File SamtoolsExe                # Path to samtools executable
   String BwaSamtoolsThreads       # Specifies the number of thread required per run
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   String DebugMode                # Flag to enable Debug Mode

   File BashPreamble               # Bash script that helps control zombie processes
   File AlignmentScript            # Bash script which is called inside the WDL script

   String AlignSoftMemLimit        # Soft memory limit - nice shutdown
   String AlignHardMemLimit        # Hard memory limit - kill immediately


   command <<<
      source ${BashPreamble}
      /bin/bash ${AlignmentScript} -s ${SampleName} -p ${Platform} -L ${Library} -f ${PlatformUnit} -c ${CenterName} -P ${PairedEnd} -l ${InputRead1} -r ${InputRead2} -G ${Ref} -e ${BWAExe} -K ${ChunkSizeInBases} -o "'${BWAExtraOptionsString}'" -S ${SamtoolsExe} -t ${BwaSamtoolsThreads} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${BwaSamtoolsThreads}"
      s_vmem: "${AlignSoftMemLimit}"
      h_vmem: "${AlignHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }

}


task mergebamTask {
   String SampleName               # Name of the Sample

   Array[File] InputBams           # Input Sorted BAM File
   Array[File] InputBais           # Input Sorted Bam Index File


   String SamtoolsExe              # Path to Samtools Executable

   String DebugMode                # Variable to check whether Debud Mode is on

   String MergeSoftMemLimit        # Soft memory limit - nice shutdown
   String MergeHardMemLimit        # Hard memory limit - kill immediately
   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File MergeBamScript             # Bash script that is called inside the WDL script

   command <<<
   	   source ${BashPreamble}
   	   /bin/bash ${MergeBamScript} -b ${sep=',' InputBams} -s ${SampleName} -S ${SamtoolsExe} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      s_vmem: "${MergeSoftMemLimit}"
      h_vmem: "${MergeHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }
}

task dedupTask {

   String SampleName               # Name of the Sample

   File InputBams                  # Input Sorted BAM File
   File InputBais                  # Input Sorted Bam Index File


   File GATKExe                    # Path to GATK4 executable
   File JavaExe                    # Path to Java8 executable
   String JavaOptionsString        # String of java vm options, like garbage collection and maximum and minimum memory. Can NOT be empty

   String DebugMode                # Variable to check whether Debud Mode is on

   String DedupSoftMemLimit        # Soft memory limit - nice shutdown
   String DedupHardMemLimit        # Hard memory limit - kill immediately
   File BashPreamble               # shell file to source before each task
   File BashSharedFunctions        # Bash script with shared functions

   File DedupScript                # Bash script that is called inside the WDL script

   command <<<
   	   source ${BashPreamble}
   	   /bin/bash ${DedupScript} -s ${SampleName} -b ${InputBams} -S ${GATKExe} -J ${JavaExe} -e ${JavaOptionsString} -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: 1
      s_vmem: "${DedupSoftMemLimit}"
      h_vmem: "${DedupHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }
}


task deliverAlignmentTask {

   File InputBams                         # aligned sorted dedupped BAM file
   File InputBais                         # aligned sorted dedupped BAM.BAI file

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File BashPreamble                      # Bash script that helps control zombie processes
   File BashSharedFunctions               # Bash script that contains shared helpful functions
   File DeliveryAlignment_Script          # Bash script that performs the delivery

   String DeliveryFolder_Alignment        # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on

   command {
      source ${BashPreamble}
      /bin/bash ${DeliveryAlignment_Script} -s ${SampleName} -b ${InputBams} -j ${WorkflowJson} -f ${DeliveryFolder_Alignment} -F ${BashSharedFunctions} ${DebugMode}
   }

}


task bqsrTask {
   String SampleName                     # Name of the Sample

   File InputBams                        # Input Sorted Deduped Bam
   File InputBais                        # Input Sorted Deduped Bam Index

   File Ref                              # Reference Genome
   File RefFai                           # Reference files that are provided as implicit inputs
   File RefDict                          # to the WDL Tool to help perform BQSR


   String BqsrKnownSites                 # List of known sites, including dbSNP
   String GenomicInterval                # Array of chromosome names or genomic intervals for parallel analysis
   File GATKExe                          # Path to GATK4 executable
   String ApplyBQSRExtraOptionsString    # String of extra options for ApplyBQSR. This can be an empty string
   File JavaExe                          # Path to Java8 executable
   String JavaOptionsString              # String of java vm options, like garbage collection and maximum and minimum memory. Can NOT be empty

   String BqsrSoftMemLimit               # Soft memory limit - nice shutdown
   String BqsrHardMemLimit               # Hard memory limit - kill immediately

   File BashPreamble                     # Bash script to run before every task
   File BashSharedFunctions              # Bash script with shared functions
   File BqsrScript                       # Path to bash script called within WDL script

   String DebugMode                      # Enable or Disable Debug Mode

   command <<<
       source ${BashPreamble}
       /bin/bash ${BqsrScript} -s ${SampleName} -b ${InputBams} -G ${Ref} -k ${BqsrKnownSites} -I ${GenomicInterval} -S ${GATKExe} -o ${ApplyBQSRExtraOptionsString} -J ${JavaExe} -e "'${JavaOptionsString}'" -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      s_vmem: "${BqsrSoftMemLimit}"
      h_vmem: "${BqsrHardMemLimit}"
   }

   output {
      File OutputBams = "${SampleName}.${GenomicInterval}.bam"
      File OutputBais = "${SampleName}.${GenomicInterval}.bai"
   }
}

task variantCallingTask {
   String SampleName                              # Name of the Sample

   File InputBams                                 # Input Sorted Deduped Bam
   File InputBais                                 # Input Sorted Deduped Bam Index

   File Ref                                       # Reference Genome
   File RefFai                                    # Reference Genome index
   File RefDict                                   # Reference Genome dictionary

   File DBSNP                                     # DBSNP file
   File DBSNPIdx                                  # Index file for the DBSNPs
   String GenomicInterval                         # Array of chromosome names or genomic intervals for parallel analysis

   File GATKExe                                   # Path to GATK4 executable
   String HaplotyperThreads
   String HaplotyperExtraOptionsString            # String of extra options for haplotyper, this can be an empty string
   File JavaExe                                   # Path to Java8 executable
   String JavaOptionsString                       # String of java vm options. Can NOT be empty

   String HaplotyperSoftMemLimit                  # Soft memory limit - nice shutdown
   String HaplotyperHardMemLimit                  # Hard memory limit - kill immediately

   File BashPreamble                              # bash script to source before every task
   File BashSharedFunctions                       # Bash script with shared functions
   File HaplotyperScript                          # Path to bash script called within WDL script

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${HaplotyperScript} -s ${SampleName} -b ${InputBams} -G ${Ref} -D ${DBSNP} -I ${GenomicInterval} -S ${GATKExe} -t ${HaplotyperThreads} -o "'${HaplotyperExtraOptionsString}'" -J ${JavaExe} -e "'${JavaOptionsString}'" -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      cpu: "${HaplotyperThreads}"
      s_vmem: "${HaplotyperSoftMemLimit}"
      h_vmem: "${HaplotyperHardMemLimit}"
   }

  output {
      File OutputVcf = "${SampleName}.${GenomicInterval}.g.vcf"
      File OutputVcfIdx = "${SampleName}.${GenomicInterval}.g.vcf.idx"
   }

}


task mergegvcfsTask {
   String SampleName                              # Name of the Sample

   Array[File] InputGvcfs                         # Input GVCF files
   Array[File] InputIdxs                          # Input GVCF files Index

   File GATKExe                                   # Path to GATK4 executable
   File JavaExe                                   # Path to Java8 executable
   String JavaOptionsString                       # String of java vm options. Can NOT be empty

   String MergeSoftMemLimit                       # Soft memory limit - nice shutdown
   String MergeHardMemLimit                       # Hard memory limit - kill immediately
   File BashPreamble                              # Bash script that helps control zombie processes
   File BashSharedFunctions                       # Bash script that contains shared helpful functions
   File MergeGvcfsScript                          # Bash script that is called inside the WDL script

   String DebugMode                               # Enable or Disable Debug Mode


   command <<<
        source ${BashPreamble}
        /bin/bash ${MergeGvcfsScript} -s ${SampleName} -b ${sep=',' InputGvcfs} -S ${GATKExe} -J ${JavaExe} -e "'${JavaOptionsString}'" -F ${BashSharedFunctions} ${DebugMode}
   >>>

   runtime {
      s_vmem: "${MergeSoftMemLimit}"
      h_vmem: "${MergeHardMemLimit}"
   }

  output {
      File OutputVcf = "${SampleName}.g.vcf"
      File OutputVcfIdx = "${SampleName}.g.vcf.idx"
   }

}

task deliverHaplotyperVCTask {

   File InputVcf                          # VCF File
   File InputVcfIdx                       # VCF.IDX File

   String SampleName                      # Name of the Sample

   File WorkflowJson                      # JSON file for the workflow

   File BashPreamble                      # Bash script that helps control zombie processes
   File BashSharedFunctions               # Bash script that contains shared helpful functions
   File DeliveryHaplotyperVC_Script       # Bash script that performs the delivery

   String DeliveryFolder_HaplotyperVC     # Path to delivery folder
   String DebugMode                       # Variable to check whether Debud Mode is on; optional

   command {
      source ${BashPreamble}
      /bin/bash ${DeliveryHaplotyperVC_Script} -s ${SampleName} -r ${InputVcf} -j ${WorkflowJson} -f ${DeliveryFolder_HaplotyperVC} -F ${BashSharedFunctions} ${DebugMode}
   }

}
