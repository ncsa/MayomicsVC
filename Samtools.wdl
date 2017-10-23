#This task calls SAM Tools. SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

task SAM_Tools {

	File InputSam
	String OutputBam	

   	command {
                 /usr/local/apps/bioapps/samtools/samtools-1.5/samtools view -@ 17 -o ${OutputBam} ${InputSam}
        }
        output {
                File out = stdout()
        }
}

workflow SAM_Tools_Run {
        call SAM_Tools
}
