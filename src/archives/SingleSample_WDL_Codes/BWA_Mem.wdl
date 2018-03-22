#This task calls BWA Mem Tool. There's a reference file and a couple of read files are the inputs
#In this step the reads are aligned with the reference data. The information is copied to a bam file.
 
task BWA_Mem {
        File RefFasta
        File Input_Read1
        File Input_Read2
        String BWA_REF
        String sampleName

        command {
                /usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} > /projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/archives/SingleSample_WDL_Codes/${sampleName}.aligned.sam
        }
        output {
                File out = "/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/src/archives/SingleSample_WDL_Codes/${sampleName}.aligned.sam"
        }

}

workflow BWA_Mem_Run {
        call BWA_Mem
}

