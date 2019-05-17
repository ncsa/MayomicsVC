# Purpose

The purpose of this readme is to understand the peculiarities of the WDL scripts. The WDL scripting structure is standard across all the scripts.

# Introduction to WDL

The Workflow Description Language (WDL) is a way to specify data processing workflows with a human-readable and -writeable syntax. WDL makes it straightforward to define analysis tasks, chain them together in workflows, and parallelize their execution. The language makes common patterns simple to express, while also admitting uncommon or complicated behavior; and strives to achieve portability not only across execution platforms, but also different types of users. Whether you are an analyst, a programmer, an operator of a production system, or any other sort of user, WDL should be accessible and understandable.

# WDL Scripting Structure

<details>
<summary>
Header
</summary>
      
```bash scripting
###########################################################################################
##              This WDL script performs alignment using BWA Mem                         ##
###########################################################################################
```

The header clearly states what action the WDL task is supposed to perform
</details>

<details>
<summary>
Variables
</summary>

```bash scripting
   File InputRead1                 # Input Read File           
   String InputRead2               # Input Read File           
   String SampleName               # Name of the Sample
      ...
   Boolean PairedEnd               # Variable to check if single ended or not
```
File InputRead1 and String InputRead2 are mentioned because they sometimes have paired ended fastqs and single ended fastqs. However, we default to single ended fastq’s. That is why the left reads are always present and Cromwell checks every variable that is a file and it will check whether that file exists or not. InputRead2 is stored as a string because we don’t want Cromwell to check it at all, in case there is no InputRead2. 
The data can be processed in three different ways, controlled by the PairedEnd variable and format of the InputRead2 variable; there are three cases:

(1) Process paired ended data as paired ended: PairedEnd=true, InputRead1 variable contains left reads, InputRead2 variable contains right reads; assume they are in correct order of lanes and files that belong to the same pair are at the same position in each list.

(2) Process paired ended data as single ended: PairedEnd=false, InputRead1 variable contains left and right reads, InputRead2 variable is empty

(3) Process single ended data as single ended: PairedEnd=false, InputRead1 variable contains left reads, InputRead2 variable is empty

```
   String Sentieon                 # Path to Sentieon
   String SentieonThreads          # Specifies the number of thread required per run

   File BashPreamble               # Bash script that helps control zombie processes
   File BashSharedFunctions        # Bash script that contains shared helpful functions
   File AlignmentScript            # Bash script which is called inside the WDL script
   File AlignEnvProfile            # File containing the environmental profile variables
   String ChunkSizeInBases         # The -K option for BWA MEM
   String BWAExtraOptionsString    # String of extra options for BWA. This can be an empty string.

   String AlignSoftMemLimit        # Soft memory limit - nice shutdown
   String AlignHardMemLimit        # Hard memory limit - kill immediately

   String DebugMode                # Flag to enable Debug Mode
 ```
</details>      

<details>
<summary>
Command Block
</summary>
      
```bash scripting
command <<<
      source ${BashPreamble}
      /bin/bash ${AlignmentScript} -P 
      ${PairedEnd} -l ${InputRead1} -r
      ${InputRead2} -s ${SampleName} -p 
      ${Platform} -L ${Library} -f 
      ${PlatformUnit} -c ${CenterName} -G 
      ${Ref} -o "'${BWAExtraOptionsString}'" -K
      ${ChunkSizeInBases} -S ${Sentieon} -t 
      ${SentieonThreads} -e ${AlignEnvProfile} 
      -F ${BashSharedFunctions} ${DebugMode}
   >>>
```

1. Bash is linked to WDL through the command block
2. The command block lists the input options with its corresponding variables and calls the shell script.
3. WDL reads the values from the json files and passes those values through the variable names defined at the top of the script into the 
shell script.
</details>

<details>
<summary>
Runtime Block      
</summary>
      
``` bash scripting
runtime {
      cpu: "${SentieonThreads}"
      s_vmem: "${AlignSoftMemLimit}"
      h_vmem: "${AlignHardMemLimit}"
   }
```

1. We need to define the soft memory limit and the hard memory limit for every task because we have one task to one bash script and one individual bash script to one individual automatic bio informatics analysis.
2. The reason why we choose to have one bioinformatics analysis per shell script and one shell script per WDL task is to avoid the complication to use different number of threads for different lines of bash scripts. 
</details>

<details>
<summary>
Output Block
</summary>
      
``` bash scripting
output {
      File OutputBams = "${SampleName}.bam"
      File OutputBais = "${SampleName}.bam.bai"
   }

} 
```
We need to have the output block because the output of one task serves as the input to other
</details>
