# The task below calls BWA Mem Tool. There's a reference file and a couple of read files are the inputs
# In this step the reads are aligned with the reference data. The information is copied to a bam file.
# In this example the BWA Mem tool is run for multiple samples

task BWA_Mem {				

        File Ref_Amb_File         
        File Ref_Dict_File
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_fai_File
        File Ref_pac_File
        File Ref_sa_File 
 
# The files declared above are additional reference files provided to the BWA Mem Tool apart from the reference Fasta File and Input Read File
	
	File RefFasta			# Reference Fasta File  
	File Input_Read1		# Input Read File 
        File Input_Read2		# Second Input Read File (Optional)
        String sampleName		# Name of the Sample

        command {
                /usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} > ${sampleName}.aligned.sam
        }
        output {
                Array[File] Aligned_Sam = glob("${sampleName}.aligned.sam")     # Output is an Array of aligned 										SAM Files
        }

}

# The task below calls SAM Tools. SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.


task Samtools {

	Array[File] aligned_sam		# Array of input aligned SAM Files
	String aligned_bam		# Output aligned BAM File

	command {
                 /usr/local/apps/bioapps/samtools/samtools-1.5/samtools view -@ 17 -o ${aligned_bam} ${sep=',' aligned_sam}
        }
        output {
                Array[File] Aligned_Bam = glob("${aligned_bam}")                # Output is an Array of aligned 									 	BAM Files
        }
}

workflow MultipleSample_MultiStepWorkFlow {
	File InputSamplesFile
        Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)
	
	scatter(sample in inputsamples)
	{
        call BWA_Mem
	{
	input :
	sampleName = sample[0],
	Input_Read1 = sample[1],
	Input_Read2 = sample[2]
	}

	call Samtools
	{
	input :
	aligned_sam = BWA_Mem.Aligned_Sam
	}
	} 
}

