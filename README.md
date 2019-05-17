# Acknowledgements

The <a href = "https://individualizedmedicineblog.mayoclinic.org/2017/08/31/the-grand-challenge-using-supercomputers-to-speed-diagnosis-and-treatment/"> Mayo Grand Challenge Project </a>  work was a product of the Mayo Clinic and Illinois Alliance for Technology-Based Healthcare. Special thanks for the funding provided by the Mayo Clinic Center for Individualized Medicine and the Todd and Karen Wanek Program for Hypoplastic Left Heart Syndrome. We also thank the Interdisciplinary Health Sciences Institute, the Carl R. Woese Institute for Genomic Biology and the National Center for Supercomputing Applications for their generous support and access to resources. We particularly acknowledge the support of Keith Stewart, M.B., Ch.B., Mayo Clinic/Illinois Grand Challenge Sponsor and Director of the Mayo Clinic Center for Individualized Medicine. Many thanks to the Sentieon team for consultation and advice on the Sentieon variant calling software.

# Steps and tools to run the workflow
### Steps
<details>
  <summary> 
   Clone the repository and load the necessary modules
  </summary>
  <br>
  a. Visit the MayomicsVC Repository and clone the repository as in the example given below:
  
  ```bash scripting
    git clone -b dev https://github.com/ncsa/MayomicsVC.git
  ```
  b. Load the necessary modules
  
This workflow requires the following tools to be installed in order to run the same, (1) Cromwell-34, (2) Java-1.8, and (3) Python-3.6.1. After installation, set the respective executable tool to the environmental $PATH

 </details>
 
 <details>
 <summary>
 Create configuration and environment profile files
  <br>
  </summary>
  <br>
  a. The user needs to provide certain input configuration files to describe the location of the data, tools, and the memory requirements, to be used in the workflow.
  
  ```bash scripting
  ## make a config directory
  mkdir Config
  cd Config
  
  ## input parameters
  touch run_info.txt
  
  ## file/software paths and options
  touch tool_info.txt
  
  ## sample names
  touch sample_info.txt
    
  ## memory info
  touch memory_info.txt
  ```
  
b. Sentieon requires a license to run. This license is a bash environmental variable, since the Sentieon commands are bash commands executed from within the pipeline. An "environmental" profile file is passed in with each task in the workflow, containing the Sentieon license environmental variable. Below are the environmental profile files that need to be present in the Config directory:

```bash scripting
ls Config/ | grep Profile

AlignEnvProfile.file 
BqsrEnvProfile.file 
DedupEnvProfile.file 
HaplotyperEnvProfile.file 
RealignEnvProfile.file
TrimEnvProfile.file 
VqsrEnvProfile.file 
```

</details>

<details>
<summary>
 Use WOM tool to create the JSON and run parser to populate the JSON
<br>
</summary>
  
a. WDL will use a json file to read in the locations data. The user first generates a json with the necessary input keys. The values will be added later.

```
mkdir Jsons
cd MayomicsVC
java -jar ${WOMTOOL} inputs MayomicsVC/src/wdl/GermlineMasterWorkflow.wdl
```
b. The JSON needs to be filled in with certain inputs and below given is an instance of the germline master workflow command:

```
cat ../Jsons/GermlineMasterWorkflow.json
{
"GermlineMasterWF.realign.RealignSoftMemLimit": "String",
"GermlineMasterWF.bqsr.DebugMode": "String",
....
"GermlineMasterWF.dedup.DedupSoftMemLimit": "String",
"GermlineMasterWF.merge.BashSharedFunctions": "File"
}
```

c. To run the parser to populate the JSON, run the following bash command:

```
 python src/python/config_parser.py -i ~/Config/run_info.txt -i ~/Config/sample_info.txt -i ~/Config/tool_info.txt --jsonTemplate ~/Jsons/<test_name>.json.tmpl -o ~/Jsons/<test_name>.json
```
</details>

<details>
<summary>
Run validator to validate entries in JSON
<br>
</summary>
Cromwell expects from the WDL file the variable types of the input variables in order to run the workflow successfully. Hence, we have written another python script to pass into the newly filled-in json file, and the key_types file from the repository:
  
```
python MayomicsVC/src/python/key_validator.py -i Jsons
/GermlineMasterWorkflow.FilledIn.json --KeyTypeFile MayomicsVC/key_types.
json
```

</details>

<details>
<summary>
Zip source code and run the script
<br>
</summary>
  
a. In order for Cromwell to know the paths to the task scripts, it is necessary to point to the scripts when executing the entire workflow and this is done by zipping the source code.

```
zip -r MayomicsVC.zip MayomicsVC
```
b. To run the germline workflow, execute the below command:  
```  
java -jar ${CROMWELL} run ./MayomicsVC/src/wdl/GermlineMasterWorkflow.wdl -i Jsons/GermlineMasterWorkflow.FilledIn.json -p MayomicsVC.zip
```
The outputs will be present in the delivery folder.
</details>

### Tools 
<details>
   <summary>
     iForge - the NCSA Industry Supercomputer <br>
  </summary>
        • Intel "Skylake" Xeon Gold 6148 <br>
        • 20-core CPU, 2.4 GHz <br>
        • Dual-CPU motherboard <br>
        • Total cores: 40 <br>
        • 192 GB RAM, 2666 MHz <br> 
        • Storage: 4+ PB <br>
        • IBM GPFS ver. 4 with custom metadata acceleration <br>
        • EDR lnfiniband, 100 GB/sec bandwidth, 100 ns latency <br>
        • WAN: 80 GB/sec <br>
        • OS: Red Hat Enterprise Linux 6 <br>
</details>

<details>
  <summary>
  BioInformatics Tools
  </summary> 
  
  <summary>
  <b> Germline</b>
  </summary>
      <a href ="https://www.sentieon.com/"> <b> Sentieon </b> </a>: sentieon-genomics-201808 software package for secondary DNA analysis
      
 <a href ="https://cutadapt.readthedocs.io/en/stable/"> <b> Cutadapt </b> </a>: software to remove adapter sequences 
  
  <summary>
  <b> Somatic </b>
  </summary>
  
  Mutect  
  Manta   
  Strelka
  
  </details>
 
 <details>
  <summary>
    Workflow Management Systems 
  </summary>
  <a href = "https://software.broadinstitute.org/wdl/"> Cromwell/WDL <a>
    </details>                                                                                             

# Objective

The objective of this project was to create the best practices genomic variant calling workflow using the workflow management system Cromwell/WDL. Our main goal was to create a workflow that is trivially maintainable in a production clinical setting or active reserach lab, so that numerous, large scale and robust analyses could be performed. Therefore, our main considerations where (1) simplicity of code, (2) robustness against misconfiguration, and (3) ease of debugging during execution.
 
# Variant Calling Workflow

The workflow has multiple components, each implemented as higher-level modules. For example, in BAM cleaning we have two modules: Alignment and Realignment/Recalibration. Each module consists of *tasks* - lowest complexity modules that represent meaningful bioinformatics processing steps (green boxes in the detailed workflow architecture below), such as alignment against a reference or deduplication of aligned BAMs. These tasks are written as .wdl scripts and are displayed below in the detailed workflow architecture:

<img src="https://user-images.githubusercontent.com/43070131/52230023-fa7b8c00-287b-11e9-82d1-2dd6146a1f3b.PNG" alt="Detailed Workflow Architecture" width="800">

<a href ="https://drive.google.com/file/d/1KpT3hou8Sb4zK4M5HzaWF2_7RLUqtWee/view?usp=sharing"> Link to edit image </a>

Description:
1. <b> Adapter Trimming </b>: Trim the adapters from the reads obtained from the sequencer using CutAdapt
2. <b> Alignment</b>: Align the reads to a reference genome using Sentieon's BWA-MEM
3. <b> Merge </b>: Merges the reads
3. <b> Mark Duplicates </b>: Remove duplicate threads
4. <b> Realignment </b>: Realign reads using Sentieon Realigner
5. <b> Base Quality Score Recalibration (BQSR) </b>: Calculate the required modification of the quality scores on the BAM produced in the Realignment stage
6. <b> Variant Caller/Haplotyper</b> : Create a VCF file containing the DNASeq variants reported for an individual
7. <b> Variant Quality Score Recalibration (VQSR) </b>: Assign a well-calibrated probability score to individual variant calls and determine the probability that sites are true

# Organization of the code

The src/ folder is broken up by language. In src/, we have 3 subfolders, (1) <a href ="https://github.com/ncsa/MayomicsVC/blob/master/src/shell/README.md"> Shell </a> - Consists of all the shell scripts that calls the bioinformatics software.  (2) <a href = "https://github.com/ncsa/MayomicsVC/blob/master/src/python/README.md"> Python </a> - Contains scripts to parse the config files for JSON and validate them for correctness. (3) <a href ="https://github.com/ncsa/MayomicsVC/blob/master/src/wdl/README.md"> WDL </a> - The shell scripts are called by WDL scripts, which perform the workflow management.
<img src="https://user-images.githubusercontent.com/43070131/57941976-f2c57d80-7895-11e9-8bae-1b6e4e9721ae.png" alt="Modularity implementation" width="800"> 

<a href ="https://drive.google.com/file/d/1KpT3hou8Sb4zK4M5HzaWF2_7RLUqtWee/view?usp=sharing">Link to edit image </a>

# Design principles

<details>
<summary>
 <b>Modularity:</b> Subdivides the workflow into individual parts independent from each other
 
 </summary>
Due to the complexity of the variant calling workflow, we break it up into modules to make it as easy to develop and maintain as possible. Thus, each bioinformatics step is its own module. WDL makes this easy by defining "tasks" and "workflows." Tasks
in our case wrap individual bioinformatics steps. These individual tasks are strung together into a master workflow: e.g. Germline or Somatic.

Below are the reasons for a modular design:
* Flexibility:
    * Can execute any part of the workflow
    * Useful for testing or after failure
    * Can swap tools in and out for every task based on user's choice
* Optimal resource utilization: can specify ideal number of nodes, walltime, etc. for every stage
* Maintainability: 
    * Can edit modules without breaking the rest of the workflow 
    * Modules like QC and user notification, which serve as plug-ins for other modules, can be changed without updating multiple places       in the workflow
The sections below explain in detail the implementation and benefits of our approach.
</details>

<details>
 <summary>
  <b> Data parallelism and scalability: </b> Parallel execution of tasks
 </summary>
Normally, the variant calling workflow must support repetitive fans and merges in the code (conditional on user choice in the runfile):
 
* Splitting of the input sequencing data into chunks, performing alignment in parallel on all chunks, 
and merging the aligned files per sample for sorting and deduplication
* Splitting of aligned/deduped BAMs for parallel realignment and recalibration per chromosome.

This is because GATK3 was not fast enough to work on a whole human genome without chunking, whereas the Sentieon variant calling implementation is very fast. Thus, we chose to keep the workflow very simple for maintainability. We do not chunk the input fastq. The workflow is implemented on a per sample basis and trimming and alignment are performed in parallel on multiple lanes. Cromwell takes care of parallelization and scalability behind the scences. We provision user control of threading and memory options for every step.
</details>

<details>
 <summary>
  <b> Real-time logging, monitoring, data provenance tracking </b>: Real time logging/monitoring progress of jobs in workflow
 </summary>

At any moment during the run, the analyst should be able to assess:

* Which stage of the workflow is running for every sample batch 
* Which samples may have failed and why 

Additionally, a well-structured post-analysis record of all events 
executed on each sample is necessary to ensure reproducibility of 
the analysis. 

Cromwell provides for most of these via the output folder structure and logs. We have added an extra layer of logging and error reporting, described below in implementation.
</details>

<details>
 <summary>
  <b> Fault tolerance and error handling </b> : Workflow should be robust against hardware/software/data failure
 
 </summary>
 
The workflow should:

* Give the user the option to fail or continue the whole workflow when something goes wrong with one of the samples
* Have the ability to move a task to a spare node in the event of hardware failure.

The latter is a function of Cromwell, but the workflow should support it by requesting a few extra nodes (beyond the nodes required based on user specifications).

To prevent avoidable failures and resource waste, the workflow should: 

* Check that all the sample files exist and have nonzero size before the workflow runs
* Check that all executables exist and have the right permissions before the workflow runs
* After running each module, check that output was actualy produced and has nonzero size
* Perform QC on each output file, write results into log, give user option to continue even if QC failed.

</details>

<details>
 <summary>
  <b> Portability </b> : Write the workflow once, deploy it in many environments.
 
 </summary>

For a workflow as complex as genomic variant calling, having to change and adapt for each different cluster is extremely 
counterproductive. Hence, the workflow should be able to port smoothly among the following three kinds of systems:
 
* grid clusters with PBS Torque
* grid clusters with OGE
* AWS

</details>

<details>
 <summary>
  <b>Development and test automation </b>: Support multiple levels of automated testing

</summary>
 
The workflow should be constructed in such a way as to support the below testing activities:
 
* Individual task testing on each task
* Integration testing for each codepath in each workflow stage
* Integration testing for the main (i.e. most used) codepath in the workflow
 </details>

# Relevant Research Papers & Sources
1. <a href = "https://wiki.ncsa.illinois.edu/display/LH/HPC+for+Computational+Genomics"> NCSA Genomics Team</a>
2. <a href = "https://www.biorxiv.org/content/10.1101/396325v1"> Computational performance and accuracy of Sentieon DNASeq variant calling workflow </a>
3. <a href = "https://www.biorxiv.org/content/10.1101/348565v1"> Performance benchmarking of GATK3.8 and GATK4 </a>
