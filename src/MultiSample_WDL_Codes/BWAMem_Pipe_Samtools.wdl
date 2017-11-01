# Pipe BWA Mem output to Samtools to create BAM File
task BWAMem_Pipe_Samtools {

        File Ref_Amb_File
        File Ref_Dict_File
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_fai_File
        File Ref_pac_File
        File Ref_sa_File
        File RefFasta                     
        File Input_Read1                 
        File Input_Read2                
        String sampleName
	String OutputBam

        command {
                $BWA bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} | /usr/local/apps/bioapps/samtools/samtools-1.5/samtools view -@ 2 -o ${OutputBam}
        }
        output {
               Array[File] Aligned_Bam = glob("${OutputBam}")   
        }

}

workflow BWAMem_Pipe_Samtools_Run {
        File InputSamplesFile
        Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)

        scatter(sample in inputsamples)
        {
        call BWAMem_Pipe_Samtools
        {
        input :
        sampleName = sample[0],
        Input_Read1 = sample[1],
        Input_Read2 = sample[2]
        }
}
}
