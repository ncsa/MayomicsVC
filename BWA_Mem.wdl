#This task calls BWA Mem Tool. There's a reference file and a couple of read files are the inputs
#In this step the reads are aligned with the reference data. The information is copied to a bam file.
 
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
        String sampleName

        command {
                /usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} > /projects/mgc/Project_1/ram/CromwellWDL/BWA_Mem/${sampleName}.aligned.sam
        }
        output {
                File out = "/projects/mgc/Project_1/ram/CromwellWDL/BWA_Mem/${sampleName}.aligned.sam"
        }

}

workflow BWA_Mem_Run {
        call BWA_Mem
}

