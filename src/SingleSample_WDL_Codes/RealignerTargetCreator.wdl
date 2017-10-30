#This task calls the Realign Target Creator  application.

 
task RealignTargetCreator {
        File RefFasta
	File FastIndex
        File DeduppedBamInput
	File GoldenVCF
	File Ref_Ann_File
	File Ref_bwt_File
	File Ref_Dict_File
	File Ref_sa_File
	File Ref_Amb_File
	File Ref_pac_File
	File DeduppedBamInputIndex
	String RealignTargetCreator_Output 

        command {
                /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar $GATK -R ${RefFasta} -I ${DeduppedBamInput} -T RealignerTargetCreator -nt 18 -known ${GoldenVCF} -o ${RealignTargetCreator_Output}
        }
        output {
                File out = stdout()
        }
}

workflow RealignTargetCreator_Run {
        call RealignTargetCreator
}

