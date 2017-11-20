# The task below calls BWA Mem Tool. There's a reference file and a couple of read files are the inputs
# In this step the reads are aligned with the reference data. The information is copied to a bam file.
# In this example the BWA Mem tool is run for multiple samples

task BWA_Mem {				

   File Ref_Amb_File         
   File Ref_Dict_File
   File Ref_Ann_File
   File Ref_bwt_File
   File Ref_fai_File
   File Ref_pac_File
   File Ref_sa_File 
 
# The files declared above are additional reference files provided to the BWA Mem Tool apart from the reference Fasta File and Input Read File
	
   File RefFasta			# Reference Fasta File  
   File Input_Read1		# Input Read File 
   File Input_Read2		# Second Input Read File (Optional)
   String sampleName		# Name of the Sample
   String Capture_Exit_Code

   command {
              /usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa mem -t 12 -M -k 32 -I 300,30 -R "@RG\tID:lane1\tLB:${sampleName}\tPL:illumina\tPU:lane1\tSM:lane1\tCN:${sampleName}" ${RefFasta} ${Input_Read1} ${Input_Read2} > ${sampleName}.aligned.sam
           
              if [ $? -ne 0 ]; then
                 echo '${sampleName} has failed at the BWA Mem Step' >> ${Capture_Exit_Code}
              fi
   }

   output {
           # Output is an Array of aligned SAM Files
              Array[File] Aligned_Sam = glob("${sampleName}.aligned.sam")     
   }

   runtime {
  	      continueOnReturnCode: true
   }
}

# The task below calls SAM Tools. SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.


task Samtools {

   Array[File] aligned_sam		# Array of input aligned SAM Files
   String sampleName		# Name of the Sample
   String Capture_Exit_Code
	
   command {
              /usr/local/apps/bioapps/samtools/samtools-1.5/samtools view -@ 17 -o ${sampleName}.aligned.bam ${sep=',' aligned_sam}
 
              if [ $? -ne 0 ]; then
                 echo '${sampleName} has failed at the Samtools Step' >> ${Capture_Exit_Code}
              fi
   }
  
   output {
           # Output is an Array of aligned BAM Files
             Array[File] Aligned_Bam = glob("${sampleName}.aligned.bam")                
   }

   runtime {
              continueOnReturnCode: true
   }
}

#This task calls the Novosort application. This tool si responsible for the sorting of the bam file. 
###Reduced run times from multi-threading and by combining sort & merge in one step.

task Novosort {
   Array[File] Aligned_Bam
   String sampleName
   String Capture_Exit_Code

   command {
              /projects/bioinformatics/builds/novocraftV3.04.06.Linux2.6/novosort -c 18 -i -o ${sampleName}.aligned.sorted.bam ${sep=',' Aligned_Bam} 
                
	      if [ $? -ne 0 ]; then
                 echo '${sampleName} has failed at the Novosort Step' >> ${Capture_Exit_Code}
              fi
   }
   
   output {
             Array[File] Aligned_Sorted_Bam = glob("${sampleName}.aligned.sorted.bam")
   }
   
   runtime {
              continueOnReturnCode: true
   }
} 

#This task calls the Picard application. This tool is responsible for the de-duplication fo the sorted bam file. 
#This tool helps in marking duplicates in the input sorted bam file

task Picard_MarkDuplicates {
   Array[File] Aligned_Sorted_Bam
   String sampleName
   String Capture_Exit_Code

   command {
              /usr/local/apps/bioapps/java/java-1.8-64bit/jre/bin/java -Xmx8g -jar /usr/local/apps/bioapps/picard/picard-2.7.1/picard.jar MarkDuplicates I=${sep=',' Aligned_Sorted_Bam } O=${sampleName}.aligned.sorted.dedupped.bam M=${sampleName}.PicardMetrics ASSUME_SORTED=true CREATE_INDEX=true

              if [ $? -ne 0 ]; then
                 echo '${sampleName} has failed at the Mark Duplicates Step' >> ${Capture_Exit_Code}
              fi
   }
   output {
             Array[File] Aligned_Sorted_Dedupped_Bam = glob("${sampleName}.aligned.sorted.dedupped.bam")
	     Array[File] PicardMetrics = glob("${sampleName}.PicardMetrics") 
   }
}

workflow AlignmentBlock {
   File InputSamplesFile
   String Capture_Exit_Code
   Array[Array[File]] inputsamples = read_tsv(InputSamplesFile)
	
   scatter(sample in inputsamples) {
      call BWA_Mem {
         input :
            sampleName = sample[0],
            Input_Read1 = sample[1],
            Input_Read2 = sample[2],
            Capture_Exit_Code = Capture_Exit_Code
      }

      call Samtools {
         input :
            sampleName = sample[0],
            aligned_sam = BWA_Mem.Aligned_Sam,
            Capture_Exit_Code = Capture_Exit_Code
      }

      call Novosort {
         input :
            sampleName = sample[0],
            Aligned_Bam = Samtools.Aligned_Bam,
            Capture_Exit_Code = Capture_Exit_Code
      }
	
      call Picard_MarkDuplicates {
         input :
            sampleName = sample[0],
            Aligned_Sorted_Bam = Novosort.Aligned_Sorted_Bam,
            Capture_Exit_Code = Capture_Exit_Code
      }

   } 
}

