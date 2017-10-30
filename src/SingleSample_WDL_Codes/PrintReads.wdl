#This program runs the Print Reads Function.

 
task Print_Reads {
        File RefFasta
	File FastIndex
	File Ref_Ann_File
	File Ref_bwt_File
	File Ref_Dict_File
	File Ref_sa_File
	File Ref_Amb_File
	File Ref_pac_File
	File Aligned_Sorted_Dedupped_Realigned_Bam
	File Aligned_Sorted_Dedupped_Realigned_BamIndex
	File BaseRecalibrator_Input
	String PrintReads_Output 

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx16g -jar $GATK -R ${RefFasta} -I ${Aligned_Sorted_Dedupped_Realigned_Bam} -T PrintReads -BQSR ${BaseRecalibrator_Input} --out ${PrintReads_Output} -nct 6
        }
        output {
                File out = stdout()
        }
} 

workflow PrintReads_Run {
        call Print_Reads
}

