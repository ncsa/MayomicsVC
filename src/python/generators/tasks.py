#!/usr/bin/env python3

from .file_types import *


class Task:
    def __init__(self, inputs, outputs, import_location, alias, task_name, delivery_task=False):
        """
        :param inputs: list of input file types
        :param outputs: list of output file types
        :param import_location: location of the task to be imported
        :param alias: the task name (in caps in the import statement, all lowercase in the task output definition)
        :param task_name: the method call from the imported task object
        :param delivery_task: if true, this task will have no child tasks; its only purpose is to save output files
        """
        self.inputs: List[FileType] = inputs
        self.outputs: List[FileType] = outputs
        self.import_location: str = import_location
        self.alias: str = alias,
        self.task_name: str = task_name,
        self.delivery_task = delivery_task,
        self.formatted_input_strings: List[str] = []
        self.dependencies_found = False

    def add_formatted_input_string(self, string):
        self.formatted_input_strings.append(string)

    def mark_dependencies_found(self):
        self.dependencies_found = True


ALIGNMENT = Task(inputs=[FASTQ],
                 outputs=[BAM, BAI],
                 import_location="src/wdl_scripts/Alignment/TestTasks/Runalignment.wdl",
                 alias="alignment",
                 task_name="RunAlignmentTask"
                 )

DEDUP = Task(inputs=[BAM, BAI],
             outputs=[BAM, BAI],
             import_location="src/wdl_scripts/Alignment/Tasks/dedup.wdl",
             alias="dedup",
             task_name="dedupTask"
             )

REALIGNMENT = Task(inputs=[BAM, BAI],
                   outputs=[BAM, BAI],
                   import_location="src/wdl_scripts/Alignment/Tasks/dedup.wdl",
                   alias="realign",
                   task_name="realignmentTask"
                   )
