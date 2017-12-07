#This task calls the Indel Realigner  application.

 
task IndelRealigner {
        File RefFasta
	File FastIndex
        File Aligned_Sorted_Dedupped_Bam
	File GoldenVCF
	File Ref_Ann_File
	File Ref_bwt_File
	File Ref_Dict_File
	File Ref_sa_File
	File Ref_Amb_File
	File Ref_pac_File
	File Aligned_Sorted_Dedupped_BamIndex
	File RealignTargetCreator_Intervals
	String IndelRealigner_Output 

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar $GATK -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Bam} -T IndelRealigner -known ${GoldenVCF} --targetIntervals ${RealignTargetCreator_Intervals} -o ${IndelRealigner_Output}
        }
        output {
                File out = stdout()
        }
}

workflow IndelRealigner_Run {
        call IndelRealigner
}

