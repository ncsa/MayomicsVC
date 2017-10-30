#This task calls GATK's tool, HaplotypeCaller in normal mode. This tool takes a pre-processed 
#bam file and discovers variant sites. These raw variant calls are then written to a vcf file.
task Haplotype_Caller {
	
	File RefFasta
        File FastIndex
        File Ref_Ann_File
        File Ref_bwt_File
        File Ref_Dict_File
        File Ref_sa_File
        File Ref_Amb_File
        File Ref_pac_File
	File GoldenVCF
	File PrintReads_Bam
        File PrintReads_BamIndex
	String HaplotypeCaller_Output

	command {
		/usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx16g -jar $GATK -T HaplotypeCaller -R ${RefFasta} --dbsnp ${GoldenVCF} -I ${PrintReads_Bam} --emitRefConfidence GVCF -gt_mode DISCOVERY --sample_ploidy 2 -nt 1 -nct 17 -o ${HaplotypeCaller_Output}

	}
	output {
		File out = stdout()
	}
}

workflow HaplotypeCaller_Run {
	call Haplotype_Caller
}
