# Background

We want a way to automatically generate Cromwell/WDL master scripts based on single words:

Example:
```
TRIMSEQ
ALIGNMENT
DELIVERY_ALIGNMENT
DEDUP
```

We created a workflow generator that takes a file with a list of tasks, and creates a workflow
  stringing inputs to the outputs of upstream tasks.
  
It works by traversing the list of tasks in reverse order, and linking the inputs of the current
  task to the closest upstream output that is of the required type (fastq, bam, vcf, etc.)
  
# Details

## Delivery tasks

Delivery tasks are nodes in the workflow that are terminating, i.e. their outputs 
  are not used further downstream. They are called purely to generate and save their
  output files

## Defining a task

The necessary information needed for each task definition can be found in tasks.py

To add new tasks, create a variable with the task's name and make an instance of the
  Task class. Additionally, add a dictionary entry in the task_dict object in tasks.py
  
An example:
```
ALIGNMENT = Task(inputs=[FASTQ],
                 outputs=[BAM, BAI],
                 import_location="src/wdl_scripts/Alignment/TestTasks/Runalignment.wdl",
                 alias="alignment",
                 task_name="RunAlignmentTask",
                 rank=2,
                 required=True
                 )
```
#### Inputs/Outputs

These are instances of the FileType class that can be found in file_types.py

Example:

`FASTQ = FileType(name="fastq", input_variable_name="InputReads", output_variable_name="OutputReads")
`

Each task lists file types in their input and output lists

#### Rank

The Rank states the order that the tasks should be executed in.
If there are multiple tools that perform the same function 
  (samtools sort, samblaster, novosort, etc.) give them the same rank.

#### Required

Some tasks in the workflow cannot be skipped, such as alignment. In this case, it must be
  defined to be a required task.
  
# Interactive Mode

Pass in the `--interactive` flag when calling the main `workflow_generator.py` script to 
  create the list of tasks in the workflow interactively.
  
Example:

```
/Users/jacobheldenbrand/Example.wdl --interactive --name ExampleWorkFlow
Select the initial task: [TRIMSEQ | ALIGNMENT | DELIVER_ALIGNMENT | DEDUP | BQSR | HAPLOTYPER]
>>> trimseq
Select the next task: [ALIGNMENT | END]
>>> alignment
Select the next task: [DELIVER_ALIGNMENT | DEDUP | BQSR | END]
>>> deliver_alignment
Select the next task: [DEDUP | BQSR | END]
>>> dedup
Select the next task: [BQSR | END]
>>> bqsr
Select the next task: [HAPLOTYPER | END]
>>> haplotyper
Select the next task: [VQSR | END]
>>> vqsr

Attempting to build a workflow with the following inputs
TRIMSEQ:		INPUTS[fastq] --> OUTPUTS[fastq]
ALIGNMENT:		INPUTS[fastq] --> OUTPUTS[bam, bai]
DELIVER_ALIGNMENT:		INPUTS[bam, bai] --> OUTPUTS[]
DEDUP:		INPUTS[bam, bai] --> OUTPUTS[bam, bai]
BQSR:		INPUTS[bam, bai] --> OUTPUTS[recal_table]
HAPLOTYPER:		INPUTS[bam, bai, recal_table] --> OUTPUTS[vcf, vcfidx]
VQSR:		INPUTS[vcf, vcfidx] --> OUTPUTS[vcf, vcfidx]

Process finished with exit code 0
```

Note, that the workflow generator does not check whether an output has been used before, it is 
  possible that the same output will be used by multiple tasks. It is up to the end user to look over
  the workflow produced and be sure that it makes sense.