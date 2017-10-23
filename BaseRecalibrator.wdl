#This task calls the Base Recalibrator  application.

 
task Base_Recalibrator {
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
	File Aligned_Sorted_Dedupped_Realigned_Bam
	File Aligned_Sorted_Dedupped_Realigned_BamIndex
	File RealignTargetCreator_Intervals
	String Base_Recalibrator_Output 

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx16g -jar $GATK -T BaseRecalibrator -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Realigned_Bam} -knownSites ${GoldenVCF} --out ${Base_Recalibrator_Output} -nct 17
        }
        output {
                File out = stdout()
        }
}

workflow Base_Recalibrator_Run {
        call Base_Recalibrator
}

