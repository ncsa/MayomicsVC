#This task calls the Novosort application. This tool si responsible for the sorting of the bam file. 
#Reduced run times from multi-threading and by combining sort & merge in one step.
 
task Novosort {
        File InputBam
        String OutputSortedBam 

        command {
                /projects/bioinformatics/builds/novocraftV3.04.06.Linux2.6/novosort -c 18 -i -o ${OutputSortedBam} ${InputBam}
        }
        output {
                File out = stdout()
        }
}

workflow Novosort_Run {
        call Novosort
}

