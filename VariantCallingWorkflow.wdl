#This wdl script is includes all the steps in the Variant Calling Workflow.

#This task calls BWA Mem Tool. There's a reference file and a couple of read files are the inputs
##In this step the reads are aligned with the reference data. The information is copied to a bam file.
 
task BWA_Mem {
        File RefFasta
        File Ref_Amb_File
        File Ref_Dict_File
	File Ref_Ann_File
	File Ref_bwt_File
	File Ref_fai_File
	File Ref_pac_File
	File Ref_sa_File
        File Input_Read1
        File Input_Read2
	String SampleName

        command {
                /usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${SampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${SampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} > ${SampleName}.sam
        }
        output {
               File aligned_sam = "${SampleName}.sam" 
        }
}

#This task calls SAM Tools. SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

task Samtools {

        File Aligned_Sam
        String Aligned_Bam

        command {
                /usr/local/apps/bioapps/samtools/samtools-1.5/samtools view -@ 17 -o ${Aligned_Bam} ${Aligned_Sam}
        }
        output {
                File aligned_bam = "${Aligned_Bam}"
        }
}

#This task calls the Novosort application. This tool si responsible for the sorting of the bam file. 
##Reduced run times from multi-threading and by combining sort & merge in one step.

task Novosort {
        File Aligned_Bam
        String aligned_sorted_bam

        command {
                /projects/bioinformatics/builds/novocraftV3.04.06.Linux2.6/novosort -c 18 -i -o ${aligned_sorted_bam} ${Aligned_Bam}
        }
        output {
                File alignedsortedbam = "${aligned_sorted_bam}"
        }
}

#This task calls the Picard application. This tool is responsible for the de-duplication fo the sorted bam file. 
##This tool helps in marking duplicates in the input sorted bam file

task Picard_MarkDuplicates {
        File aligned_sorted_bam
        String aligned_sorted_dedupped_bam
        String PicardMetrics

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar /usr/local/apps/bioapps/picard/picard-2.7.1/picard.jar MarkDuplicates I=${aligned_sorted_bam} O=${aligned_sorted_dedupped_bam} M=${PicardMetrics} ASSUME_SORTED=true CREATE_INDEX=true
        }
        output {
                File AlignedSortedDeduppedBam= "${aligned_sorted_dedupped_bam}"
        }
}

#This task calls the Realigner Target Creator  application.

task RealignerTargetCreator {
        File RefFasta
        File AlignedSortedDeduppedBam
        File GoldenVCF
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_Dict_File
        File Ref_sa_File
        File Ref_Amb_File
        File Ref_pac_File
	File Ref_fai_File
        File AlignedSortedDeduppedBamIndex
        String intervals

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar $GATK -R ${RefFasta} -I ${AlignedSortedDeduppedBam} -T RealignerTargetCreator -nt 18 -known ${GoldenVCF} -o ${intervals}
        }
        output {
                File Intervals = "${intervals}"
        }
}

#This task calls the Indel Realigner  application.

task IndelRealigner {
        File RefFasta
        File AlignedSortedDeduppedBam
        File GoldenVCF
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_Dict_File
	File Ref_fai_File
        File Ref_sa_File
        File Ref_Amb_File
        File Ref_pac_File
        File AlignedSortedDeduppedBamIndex
        File RealignTargetCreator_Intervals
        String AlignedSortedDeduppedRealignedBam

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar $GATK -R ${RefFasta} -I ${AlignedSortedDeduppedBam} -T IndelRealigner -known ${GoldenVCF} --targetIntervals ${RealignTargetCreator_Intervals} -o ${AlignedSortedDeduppedRealignedBam}
        }
        output {
                File IndelRealign_out = "${AlignedSortedDeduppedRealignedBam}"
        }
}

#This task calls the Base Recalibrator  application.

task Base_Recalibrator {
        File RefFasta
        File Ref_fai_File
        File GoldenVCF
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_Dict_File
        File Ref_sa_File
        File Ref_Amb_File
        File Ref_pac_File
        File Aligned_Sorted_Dedupped_Realigned_Bam
	File Aligned_Sorted_Dedupped_Realigned_BamIndex
        String recal_report_grp

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx16g -jar $GATK -T BaseRecalibrator -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Realigned_Bam} -knownSites ${GoldenVCF} --out ${recal_report_grp} -nct 17
        }
        output {
                File Base_Recalibrator_out = "${recal_report_grp}"
        }
}

#This program runs the Print Reads Function.

task Print_Reads {
        File RefFasta
        File Ref_fai_File
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_Dict_File
        File Ref_sa_File
        File Ref_Amb_File
        File Ref_pac_File
        File Aligned_Sorted_Dedupped_Realigned_Bam
        File Aligned_Sorted_Dedupped_Realigned_BamIndex
        File BaseRecalibrator_Input
        String aligned_sorted_dedupped_realigned_recalibrated_bam

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx16g -jar $GATK -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Realigned_Bam} -T PrintReads -BQSR ${BaseRecalibrator_Input} --out ${aligned_sorted_dedupped_realigned_recalibrated_bam} -nct 6
        }
        output {
                File PReads_out = "${aligned_sorted_dedupped_realigned_recalibrated_bam}"
        }
}

#This task calls GATK's tool, HaplotypeCaller in normal mode. This tool takes a pre-processed 
#bam file and discovers variant sites. These raw variant calls are then written to a vcf file.

task Haplotype_Caller {

        File RefFasta
        File Ref_fai_File
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_Dict_File
        File Ref_sa_File
        File Ref_Amb_File
        File Ref_pac_File
        File GoldenVCF
        File Aligned_Sorted_Dedupped_Realigned_Recalibrated_Bam
        File Aligned_Sorted_Dedupped_Realigned_Recalibrated_BamIndex
        String HaplotypeCaller_Output

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx16g -jar $GATK -T HaplotypeCaller -R ${RefFasta} --dbsnp ${GoldenVCF} -I ${Aligned_Sorted_Dedupped_Realigned_Recalibrated_Bam} --emitRefConfidence GVCF -gt_mode DISCOVERY --sample_ploidy 2 -nt 1 -nct 17 -o ${HaplotypeCaller_Output}

        }
        output {
                File HaplotypeCaller_out = "${HaplotypeCaller_Output}"
        }
}

workflow VariantCallingWorkflow_run {

        File RefFasta
        File Ref_Amb_File
        File Ref_Dict_File
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_fai_File
        File Ref_pac_File
        File Ref_sa_File
        File Input_Read1
        File Input_Read2
        String SampleName
	File AlignedSortedDeduppedBamIndex
	File Aligned_Sorted_Dedupped_Realigned_BamIndex
	File Aligned_Sorted_Dedupped_Realigned_Recalibrated_BamIndex
	File GoldenVCF

	call BWA_Mem {
	input:
        RefFasta = RefFasta,
        Ref_Amb_File = Ref_Amb_File,
        Ref_Dict_File = Ref_Dict_File, 
        Ref_Ann_File =  Ref_Ann_File,
        Ref_bwt_File = Ref_bwt_File,
        Ref_fai_File = Ref_fai_File,
        Ref_pac_File = Ref_pac_File,
        Ref_sa_File = Ref_sa_File,	
	Input_Read1 = Input_Read1,
	Input_Read2 = Input_Read2, 
	SampleName = SampleName	
	}

	call Samtools {
	input:
	Aligned_Sam = BWA_Mem.aligned_sam
	}
	
	call Novosort {
	input:
	Aligned_Bam = Samtools.aligned_bam
	}
	
	call Picard_MarkDuplicates {
	input:
	aligned_sorted_bam = Novosort.alignedsortedbam
	}
	
	call RealignerTargetCreator {
	input:
        RefFasta = RefFasta,
        Ref_Amb_File = Ref_Amb_File,
        Ref_Dict_File = Ref_Dict_File,
        Ref_Ann_File =  Ref_Ann_File,
        Ref_bwt_File = Ref_bwt_File,
        Ref_fai_File = Ref_fai_File,
        Ref_pac_File = Ref_pac_File,
        Ref_sa_File = Ref_sa_File,
	AlignedSortedDeduppedBamIndex = AlignedSortedDeduppedBamIndex,
	AlignedSortedDeduppedBam = Picard_MarkDuplicates.AlignedSortedDeduppedBam,
	GoldenVCF = GoldenVCF
	}
	
	call IndelRealigner {
	input:
        RefFasta = RefFasta,
        Ref_Amb_File = Ref_Amb_File,
        Ref_Dict_File = Ref_Dict_File,
        Ref_Ann_File =  Ref_Ann_File,
        Ref_bwt_File = Ref_bwt_File,
        Ref_fai_File = Ref_fai_File,
        Ref_pac_File = Ref_pac_File,
        Ref_sa_File = Ref_sa_File,
        GoldenVCF = GoldenVCF,
	AlignedSortedDeduppedBamIndex = AlignedSortedDeduppedBamIndex,
	AlignedSortedDeduppedBam = Picard_MarkDuplicates.AlignedSortedDeduppedBam,
        RealignTargetCreator_Intervals = RealignerTargetCreator.Intervals
	}
	
	call Base_Recalibrator {
	input:
        RefFasta = RefFasta,
        Ref_Amb_File = Ref_Amb_File,
        Ref_Dict_File = Ref_Dict_File,
        Ref_Ann_File =  Ref_Ann_File,
        Ref_bwt_File = Ref_bwt_File,
        Ref_fai_File = Ref_fai_File,
        Ref_pac_File = Ref_pac_File,
        Ref_sa_File = Ref_sa_File,
        GoldenVCF = GoldenVCF,
        Aligned_Sorted_Dedupped_Realigned_Bam = IndelRealigner.IndelRealign_out,
        Aligned_Sorted_Dedupped_Realigned_BamIndex = Aligned_Sorted_Dedupped_Realigned_BamIndex
	}
	
	call Print_Reads {
	input:
        RefFasta = RefFasta,
        Ref_Amb_File = Ref_Amb_File,
        Ref_Dict_File = Ref_Dict_File,
        Ref_Ann_File =  Ref_Ann_File,
        Ref_bwt_File = Ref_bwt_File,
        Ref_fai_File = Ref_fai_File,
        Ref_pac_File = Ref_pac_File,
        Ref_sa_File = Ref_sa_File,
        Aligned_Sorted_Dedupped_Realigned_Bam = IndelRealigner.IndelRealign_out,
        Aligned_Sorted_Dedupped_Realigned_BamIndex = Aligned_Sorted_Dedupped_Realigned_BamIndex
	}

	call Haplotype_Caller {
	input:
        RefFasta = RefFasta,
        Ref_Amb_File = Ref_Amb_File,
        Ref_Dict_File = Ref_Dict_File,
        Ref_Ann_File =  Ref_Ann_File,
        Ref_bwt_File = Ref_bwt_File,
        Ref_fai_File = Ref_fai_File,
        Ref_pac_File = Ref_pac_File,
        Ref_sa_File = Ref_sa_File,
        GoldenVCF = GoldenVCF,
	Aligned_Sorted_Dedupped_Realigned_Recalibrated_Bam = Print_Reads.PReads_out,
        Aligned_Sorted_Dedupped_Realigned_Recalibrated_BamIndex = Aligned_Sorted_Dedupped_Realigned_Recalibrated_BamIndex
	}
	 
}
