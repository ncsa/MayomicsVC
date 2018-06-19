###########################################################################################

##                     This WDL script performs BQSR using Sentieon                     ##

##                                Script Options
#       -t        "Number of Threads"                                     (Optional)
#       -G        "Reference Genome"                                      (Required)
#       -b        "Input Deduped Bam"                                     (Required)
#       -k        "List of Known Sites"                                   (Required)
#       -D        "Path to the DBSNP File"                                (Required)
#       -s        "Name of the sample"                                    (Optional)
#       -S        "Path to the Sentieon Tool"                             (Required)
#       -O        "Directory for the Output"                              (Required)
#       -L        "Sentieon License File"                                 (Required)
#       -e        "Path to the Error Log File"                            (Required)
#       -d        "Debug Mode Toggle"                                     (Optional)

############################################################################################

task bqsrTask {

   File InputBam                   # Input Sorted Deduped Bam
   File InputBamIdx                # Input Sorted Deduped Bam Index
   File RefFasta                   # Reference Genome

   File Ref_Amb_File               # These are reference files that are provided as implicit inputs
   File Ref_Fai_File               # to the WDL Tool to help perform the realignment

   String sampleName               # Name of the Sample

   File Known_Sites                # List of known sites
   File KnownSitesIdx              # Index file for the known sites

   File DBSNP                      # DBSNP file
   File DBSNP_Idx                  # Index file for the DBSNPs   

   String Threads                  # No of Threads for the Tool
   String Sentieon_License         # Sentieon License Information
   String Sentieon                 # Path to Sentieon
   String OutDir                   # Output Directory
   Boolean Debug_Mode_EN           # Enable or Disable Debug Mode
   String Error_Logs               # Filepath to Error Logs

   File BQSR_Script                # Path to bash script called within WDL script

   command {

      /bin/bash ${BQSR_Script} -s ${sampleName} -O ${OutDir} -S ${Sentieon} -r ${RefFasta} -t ${Threads} -b ${InputBam} -D ${DBSNP} -k ${Known_Sites} -e ${Error_Logs} -d ${Debug_mode_EN}

   }

   
   output {

      File RecalTable = "${OutDir}/${sampleName}.recal_data.table"
      File RecalTablePost = "${OutDir}/${sampleName}..recal_data.table.post"
      File RecalCSV = "${OutDir}/${sampleName}.recal.csv"
      File RecalPlots = "${OutDir}/${sampleName}.recal_plots.pdf"

   }

}
