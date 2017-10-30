#This task calls the Picard application. This tool is responsible for the de-duplication fo the sorted bam file. 
#This tool helps in marking duplicates in the input sorted bam file
 
task Picard_MarkDuplicates {
        File InputSortedBam
        String DeduppedOutputBam
	String PicardMetrics 

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar /usr/local/apps/bioapps/picard/picard-2.7.1/picard.jar MarkDuplicates I=${InputSortedBam} O=${DeduppedOutputBam} M=${PicardMetrics} ASSUME_SORTED=true CREATE_INDEX=true
        }
        output {
                File out = stdout()
        }
}

workflow Picard_Run {
        call Picard_MarkDuplicates
}

