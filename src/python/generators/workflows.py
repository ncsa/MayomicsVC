#!/usr/bin/env python3

import sys
from .tasks import *
from typing import List


class Workflow:
    def __init__(self, workflow_name, task_list_file, logger: ProjectLogger):
        self.workflow_name = workflow_name
        self.task_list: List[Task] = Workflow.__process_task_list_file(task_list_file)

    @staticmethod
    def __process_task_list_file(file_path: str):
        with open(file_path, "r") as F:
            trimmed_lines = [line.strip() for line in F]
        return trimmed_lines

    @staticmethod
    def __format_import_statement(node: Task):
        return 'import "' + node.import_location + '" as ' + node.alias.upper()


    def __format_task_statement():
        pass

    def __find_dependencies(self, task, task_index):
        """
        Searches through the task list in reverse order, searching for the inputs to the current task

        The function searches each upstream task and finds the first output that matches the required type that it needs
          (that is not marked as a "delivery" task)

        If the task is the first in the workflow, it has no workflow dependencies
        Fails if no outputs of the required type are found

        Side effect of marking dependencies as found on the Task itself, and saves the formatted input strings needed to
          write the task call to the final output file

        :param task: the current task
        :param task_index: its position in the task list
        """

        if task_index == 0:
            # This is the first task, by definition, it has no dependencies
            task.mark_dependencies_found()
        else:
            current_task = task_index - 1
            while current_task >= 0:
                pass



    def processTaskList(self):
        for index, task in enumerate(self.task_list):
            find_dependencies(task, index)



def main(args):
    a: List[str] = ("po", "ro")
    print(a.inputs)


if __name__ == "__main__":
    main(sys.argv[1:])
