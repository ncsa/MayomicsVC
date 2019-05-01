# Acknowledgements

This work was a product of the Mayo Clinic and Illinois Strategic Alliance for Technology-Based Healthcare. Special thanks for the funding provided by the Mayo Clinic Center for Individualized Medicine and the Todd and Karen Wanek Program for Hypoplastic Left Heart Syndrome. We also thank the Interdisciplinary Health Sciences Institute, UIUC Institute for Genomic Biology and the National Center for Supercomputing Applications for their generous support and access to resources. We particularly acknowledge the support of Keith Stewart, M.B., Ch.B., Mayo Clinic/Illinois Grand Challenge Sponsor and Director of the Mayo Clinic Center for Individualized Medicine. Many thanks to the Sentieon team for consultation and advice on the Sentieon variant calling software.

# Quick steps to run the workflow

# Objective

The objective of this project was to create the best practices genomic variant calling workflow using the workflow management system Cromwell/WDL. Our main goal was to create a workflow that is trivially maintainable in a production clinical setting or active reserach lab, so that numerous, large scale and robust analyses could be performed. Therefore, our main considerations where (1) simplicity of code, (2) robustness against misconfiguration, and (3) ease of debugging during execution.
 
# Variant Calling Workflow

The workflow has multiple components, each implemented as higher-level modules.For example, in BAM cleaning we have two modules: Alignment and Realignment/Recalibration.Each module consists of *tasks* - lowest complexity modules that represent meaningful bioinformatics processing steps (green boxes in the detailed workflow architecture below), such as alignment against a reference or deduplication of aligned BAMs.These tasks are written as .wdl scripts and is displayed below in the detailed workflow architecture:

<img src="https://user-images.githubusercontent.com/43070131/52230023-fa7b8c00-287b-11e9-82d1-2dd6146a1f3b.PNG" alt="Detailed Workflow Architecture" width="800">

<a href ="https://drive.google.com/file/d/1KpT3hou8Sb4zK4M5HzaWF2_7RLUqtWee/view?usp=sharing"> Link to edit image </a>

Steps:
1. Adapter Trimming: Trim the adapters from the reads obtained from the sequencer using CutAdapt
2. Alignment: Align the reads to a reference genome using Sentieon's BWA-MEM
3. Merge: Merge's the reads
3. Mark Duplicates: Remove duplicate threads
4. Realignment: Realign reads using Sentieon Realigner
5. Base Quality Score Recalibration (BQSR): Calculate the required modification of the quality scores on the BAM produced in the Realignment stage
6. Variant Caller/Haplotyper: Create a VCF file containing the DNASeq variants reported for an individual
7. Variant Quality Score Recalibration (VQSR): Assign a well-calibrated probability score to individual variant calls and determine the probability that sites are true

# Organization of the code

The src/ folder is broken up by language. In src/, we have 3 subfolders, (1) <a href ="https://github.com/ncsa/MayomicsVC/blob/master/src/shell/README.md"> Shell </a> - The folder shell consists of all the shell scripts that calls the bioinformatics software.  (2) Python - The python folder contains scripts to parse the config files for JSON and validate them for correctness. (3) WDL - The shell script is called by WDL which are the scripts of the workflow management.
<img src="https://user-images.githubusercontent.com/43070131/52072925-c8042300-254b-11e9-9ea8-42fa71aaa15e.PNG" alt="Modularity implementation" width="800"> 

<a href ="https://drive.google.com/file/d/1KpT3hou8Sb4zK4M5HzaWF2_7RLUqtWee/view?usp=sharing">Link to edit image </a>

# Tools to run the workflow


# Detailed steps to run the workflow
Below given are the steps to run the workflow:
<details>
  <summary> 
    1. Clone the repository and load the necessary modules
  </summary>
  
  a. Visit the MayomicsVC Repository and clone the repository as the sample example given below:
  
  ```bash scripting
  ## Current working directory
     pwd
    /projects/abv/mken/variant_calling_demo
  ```
  ```bash scripting
   ## Clone the repo
    ## Most up to date branch is currently dev
    git clone -b dev https://github.com/ncsa/MayomicsVC.git
  ```
  a. Load the necessary modules
  This workflow requires the cromwell execution engine and Java run
  
  ```bash scripting
  module load /usr/local/apps/bioapps/modules/java/java-1.8
  module load /usr/local/apps/bioapps/modules/cromwell/cromwell-34
  module load /usr/local/apps/bioapps/modules/python/python-3.6.1
  ```
 </details>
 
 <details>
 <summary>
2. Create configuration files
  </summary>
  
  The user needs to provide certain input configuration files to describe the location of the data, tools, and the memory requirements,   to be used in the workflow.
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
  </details>
  
  <details>
  <summary>
3. Create environmental profile files
  </summary>
  
Senteion requires a license to run. This license is a bash environmental variable, since the Senteion commands are bash commands executed from within the pipeline. An "environmental" profile file is passed in with each task in the workflow, containing the Senteion license environmental variable. The user defined the names of these files in the tool_info.txt config file. For organization purposes, these files should be in the Config directory that was created earlier. The liscense on iForge is used in this example. Following are the necessary environmental profiles in the Config dir:

ls Config/ | grep Profile

AlignEnvProfile.file
BqsrEnvProfile.file
DedupEnvProfile.file
HaplotyperEnvProfile.file
RealignEnvProfile.file
TrimEnvProfile.file
VqsrEnvProfile.file

Each file contains the same thing:

cat Config/AlignEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

cat BqsrEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

cat DedupEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

cat HaplotyperEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

cat RealignEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

cat TrimEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

cat VqsrEnvProfile.file
export SENTIEON_LICENSE=bwlm3.ncsa.illinois.edu:8989

</details>

<details>
<summary>
4. Use WOM tool to create JSON
</summary>
  
WDL will use a json file to read in the locations data. The user first generates a json with the necessary input keys. The values will be added later.

```
mkdir Jsons
cd MayomicsVC
java -jar ${WOMTOOL} inputs src/wdl/GermlineMasterWorkflow.wdl > ../Jsons
/GermlineMasterWorkflow.json
```
The JSON needs to be filled in with the below commands

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
</details>  
<details>
<summary>
5. Run parser to populate JSON
</summary>

Go to MayomicsVC and run the following bash command

```
 python src/python/config_parser.py -i ~/Config/run_info.txt -i ~/Config/sample_info.txt -i ~/Config/tool_info.txt --jsonTemplate ~/Jsons/<test_name>.json.tmpl -o ~/Jsons/<test_name>.json
```
</details>

<details>
<summary>
6. Run validator to validate entries in JSON
</summary>
In order for the workflow to run successfully, the variable types of the input variables must be what the Cromwell expects from the WDL code. We have writted another python script to ensure that this is the case.Pass in the newly filled in json file, and the key_types file from the repository:
  
```
python MayomicsVC/src/python/key_validator.py -i Jsons
/GermlineMasterWorkflow.FilledIn.json --KeyTypeFile MayomicsVC/key_types.
json
```

</details>

<details>
<summary>
 7. Zip source code </summary>
  
When calling tasks from within workflows, one has to use the "import" statement and explicitly refer to the task using the specific folder path leading to it. In order for Cromwell to know the paths of the task scripts, it is necessary to point to the scripts when executing the entire workflow.
This is done by passing in a zip archive containing all the scripts in there respective directories with the -p option when running the workflow (This archive will be used when executing the whole workflow later). Since the WDL code is written with the known locations of the task scripts in the repository, you can simply zip the files within the MayomicsVS. Make sure to cd into the directory before zipping though, or else the file pahts will not be corect.The workflow won't run correctly if the zip file is created in the wrong directory.

```
cd MayomicsVC
zip -r MayomicsVC.zip ./
mv MayomicsVC.zip ../
cd ../
```
</details>

<details>
<summary>
8. Viewing Outputs
</summary>
  
The outputs are in the delivery folders. From the Alignment Block, a BAM is produced, and from the HaplotyperVC block, a VCF and index is produced:

```
ls Delivery/Alignment/
GermlineMasterWorkflow.FilledIn.json NEAT_synthetic.bam NEAT_synthetic.
bam.bai
ls Delivery/HaplotyperVC/
GermlineMasterWorkflow.FilledIn.json NEAT_synthetic.vcf NEAT_synthetic.
vcf.idx
```
</details>

<details>
<summary>
9. Running the script </summary>
  
```  
java -jar $CROMWELL run <full_path_to_wdl_file>.wdl -i ~/Jsons/<test_name>.json -p MayomicsVC.zip
java -jar $WOMTOOL inputs src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl > ~/Jsons/TestTrimSequences.json.tmpl
java -jar $CROMWELL run MayomicsVC/src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl -i ~/Jsons/TestTrimSequences.json -p MayomicsVC.zip
```

</details>



# Design principles

<details>
<summary>
Modularity
 </summary>
Due to the complexity of the variant calling workflow, we break it up into modules to make it as easy to develop and maintain as possible.Thus, each bioinformatics step is its own module.WDL makes this easy by defining "tasks" and "workflows." Tasks
in our case wrap individual bioinformatics steps. These individual tasks are strung together into a master workflow: e.g. Germline or Somatic.

Below given are the reasons for a modular design:
* Flexibility:
    * Can execute any part of the workflow
    * Useful for testing or after failure
    * Can swap tools in and out for every task based on user's choice
* Optimal resource utilization: can specify ideal number of nodes, walltime, etc. for every stage
* Maintainability: 
    * Can edit modules without breaking the rest of the workflow 
    * Modules like QC and user notification, which serve as plug-ins for other modules, can be changed without updating multiple places       in the workflow
The sections below explain in detial the implementation and benefits of our approach.
</details>

<details>
 <summary>
 Data parallelism and scalability
 </summary>
Normally, the variant calling workflow must support repetitive fans and merges in the code (conditional on user choice in the runfile):
* Splitting of the input sequencing data into chunks, performing alignment in parallel on all chunks, 
and merging the aligned files per sample for sorting and deduplication
* Splitting of aligned/dedupped BAMs for parallel realignment and recalibration per chromosome.
This is because GATK3 was not fast enough to work on a whole human genome without chunking whereas GATK4 already runs faster without chunking the data and will be faster still in the future. Additionally, the Sentieon implementation is very fast as well. Thus, we chose to keep the workflow very simple for maintainability.We do not chunk the input fastq. The workflow is implemented on a per sample basis and trimming and alignment is performed in parallel on multiple lanes. Cromwell takes care of parallelization and scalability behind the scences. We provision user control of threading and memory options for every step.
</details>

<details>
 <summary>
 Real-time logging and monitoring, data provenance tracking
 </summary>
 
The workflow should have a good system for logging and monitoring progress of the jobs. 
At any moment during the run, the analyst should be able to assess: 
* Stage of the workflow is running for every sample batch 
* Samples may have failed and why 

Additionally, a well-structured post-analysis record of all events 
executed on each sample is necessary to ensure reproducibility of 
the analysis. 

Cromwell provides for most of these via the output folder structure and logs. We have added an extra layer of logging and error reporting, described below in implementation.
</details>

<details>
 <summary>
 Fault tolerance and error handling
 </summary>
 
The workflow must be robust against hardware/software/data failure. It should:
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
 Portability
 </summary>
The aim of this design principle is that a developer should be able to write a workflow once and then deploy it in many 
environments. For a workflow as complex as genomic variant calling, having to change and adapt for each different cluster is extremely 
counterproductive. Hence, the workflow should be able to port smoothly among the following three kinds of systems:
* grid clusters with PBS Torque
* grid clusters with OGE
* AWS
</details>

<details>
 <summary>
 Development and test automation 
</summary>
The workflow should be constructed in such a way as to support multiple levels of automated testing:
* Individual task testing on each task
* Integration testing for each codepath in each workflow stage
* Integration testing for the main (i.e. most used) codepath in the workflow
 </details>
