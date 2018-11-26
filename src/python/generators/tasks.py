#!/usr/bin/env python3

from .file_types import *
from typing import List
from typing import Dict


class Task:
    def __init__(self,
                 inputs,
                 outputs,
                 import_location,
                 alias,
                 task_name,
                 rank: int,
                 delivery_task=False,
                 required=False
                 ):
        """
        :param inputs: list of input file types
        :param outputs: list of output file types
        :param import_location: location of the task to be imported
        :param alias: the task name (in caps in the import statement, all lowercase in the task output definition)
        :param task_name: the method call from the imported task object
        :param rank: the rank of the task if they all were placed in order
        :param delivery_task: if true, this task will have no child tasks; its only purpose is to save output files
        :param required: this task must be included in all workflows
        """
        self.inputs: List[FileType] = inputs
        self.outputs: List[FileType] = outputs
        self.import_location: str = import_location
        self.alias: str = alias
        self.task_name: str = task_name
        self.formatted_input_strings: List[str] = []
        self.rank: int = rank

        self.delivery_task = delivery_task
        self.required = required

        self.dependencies_found = False

    def add_formatted_input_string(self, string):
        self.formatted_input_strings.append(string)

    def mark_dependencies_found(self):
        self.dependencies_found = True

    def get_input_output_summary_string(self):
        def format_list(io_types: List[FileType]):
            if len(io_types) == 0:
                return "[]"
            else:
                built_string = "["
                for i in io_types[:-1]:
                    built_string = built_string + i.name + ", "
                built_string = built_string + io_types[-1].name + "]"
                return built_string

        return (self.alias.upper() + ":\t\tINPUTS" + format_list(self.inputs) +
                " --> " + "OUTPUTS" + format_list(self.outputs)
                )


###
# Delivery tasks should not have outputs
###


CUTADAPTTRIM = Task(inputs=[FASTQ],
                    outputs=[FASTQ],
                    import_location="src/wdl_scripts/Alignment/TestTasks/Runtrim_sequences.wdl",
                    alias="trimseq",
                    task_name="RunTrimSequencesTask",
                    rank=1
                    )

ALIGNMENT = Task(inputs=[FASTQ],
                 outputs=[BAM, BAI],
                 import_location="src/wdl_scripts/Alignment/TestTasks/Runalignment.wdl",
                 alias="alignment",
                 task_name="RunAlignmentTask",
                 rank=2,
                 required=True
                 )

DELIVER_ALIGNMENT = Task(inputs=[BAM, BAI],
                         outputs=[],
                         import_location="src/wdl_scripts/DeliveryOfAlignment/Tasks/deliver_alignment.wdl",
                         alias="deliver_alignment",
                         task_name="deliverAlignmentTask",
                         delivery_task=True,
                         rank=3
                         )

DEDUP = Task(inputs=[BAM, BAI],
             outputs=[BAM, BAI],
             import_location="src/wdl_scripts/Alignment/Tasks/dedup.wdl",
             alias="dedup",
             task_name="dedupTask",
             rank=4
             )

REALIGNMENT = Task(inputs=[BAM, BAI],
                   outputs=[BAM, BAI],
                   import_location="src/wdl_scripts/Alignment/Tasks/dedup.wdl",
                   alias="realign",
                   task_name="realignmentTask",
                   rank=5
                   )

BQSR = Task(inputs=[BAM, BAI],
            outputs=[RECAL_TABLE],
            import_location="src/wdl_scripts/HaplotyperVC/Tasks/bqsr.wdl",
            alias="bqsr",
            task_name="bqsrTask",
            rank=6
            )

HAPLOTYPER = Task(inputs=[BAM, BAI, RECAL_TABLE],
                  outputs=[VCF, VCFIDX],
                  import_location="src/wdl_scripts/HaplotyperVC/Tasks/haplotyper.wdl",
                  alias="haplotyper",
                  task_name="variantCallingTask",
                  rank=7
                  )

VQSR = Task(inputs=[VCF, VCFIDX],
            outputs=[VCF, VCFIDX],
            import_location="src/wdl_scripts/HaplotyperVC/Tasks/vqsr.wdl",
            alias="vqsr",
            task_name="vqsrTask",
            rank=8
            )

task_dict: Dict[str, FileType] = {
    "CUTADAPTTRIM": CUTADAPTTRIM,
    "ALIGNMENT": ALIGNMENT,
    "DEDUP": DEDUP,
    "DELIVER_ALIGNMENT": DELIVER_ALIGNMENT,
    "BQSR": BQSR,
    "HAPLOTYPER": HAPLOTYPER,
    "VQSR": VQSR
}
