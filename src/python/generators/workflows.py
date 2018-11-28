#!/usr/bin/env python3

import sys
from .tasks import *
from typing import List
import logging
from util.log import ProjectLogger

"""
Exit code Rules:

1. Exit codes in this module are only given when an error has occurred, so they will all start with 'E.'
2. The letters 'wfg.' because they are coming from the workflow generator component of the code
3. A three letter code that hints at what the problem was
4. A number that can help to differentiate similar error codes


Error Code List
========================================================================================================================
E.wfg.Bad.1 = Workflow was created without passing in either a task list or a task list file
E.wfg.InO.1 = A task could not infer the location of one of its inputs from any upstream tasks
E.wfg.BIn.1 = An input task name was not a valid task
E.wfg.Ord.1 = The task list had tasks that were out of order
"""


class Workflow:

    INDENT = "   "

    def __init__(self, workflow_name, debug_mode, job_id="NA", input_file: str=None, task_list: List[Task]=None):
        self.workflow_name = workflow_name

        # Initialize the project logger
        if debug_mode:
            self.project_logger = ProjectLogger(job_id, "generators.workflows.Workflow", logging.DEBUG)
        else:
            self.project_logger = ProjectLogger(job_id, "generators.workflows.Workflow")
        self.job_id = job_id

        # Either input_file or task_list must be set
        if input_file is not None:
            self.task_list: List[Task] = Workflow.__process_task_list_file(self, input_file)
        elif task_list is not None:
            self.task_list: List[Task] = task_list
        else:
            self.project_logger.log_error("E.wfg.Bad.1", "Must use either interactive mode or pass in an input file")

    def __validate_task_order(self, tasks):
        """
        Given a list of tasks, assures that their ranks are in ascending order

        If tasks are out of order, an error is thrown

        :param tasks: A list of Task objects
        """
        for i in range(0, len(tasks) - 2):
            task_1 = tasks[i]
            task_2 = tasks[i + 1]

            if task_1.rank >= task_2.rank:
                self.project_logger.log_error(
                    "E.wfg.Ord.1",
                    "Cannot have " + task_1.alias.upper() + "before " + task_2.alias.upper()
                )
                sys.exit(1)

    def __process_task_list_file(self, file_path: str) -> List[FileType]:
        """
        Read in a file with one task name per line, and verify that the lines are valid task names, and convert the task
          name into its corresponding FileType instance

        :param file_path: The path to the file
        :return: a list tasks
        """
        with open(file_path, "r") as F:
            trimmed_lines = [line.strip() for line in F]

        tasks = []

        for i in trimmed_lines:
            upper_i = i.upper()
            if upper_i in task_dict.keys():
                tasks.append(task_dict.get(upper_i))
            else:
                self.project_logger.log_error("E.wfg.BIn.1", "The input task : " + i + " is not a valid task name")
                sys.exit(1)

        # Verify that the tasks are not out of order (throw an error otherwise)
        self.__validate_task_order(tasks)
        return tasks

    @staticmethod
    def __format_import_statement(node: Task):
        return 'import "' + node.import_location + '" as ' + node.alias.upper()

    @staticmethod
    def __format_import_string(input_type: FileType, output_type: FileType, task_containing_output: Task):
        """
        Takes in the types of the input and output files and the task that contains the output file and returns a string
          with the output assigned to the input variable

        Example:    InputBams = align.OutputBams

        :param input_type
        :param output_type
        :param task_containing_output The task that
        :return: Returns a string in the input format for a WDL task (InputBams = align.OutputBams)
        """
        return input_type.input_variable_name + " = " + task_containing_output.alias.lower() + "." + \
            output_type.output_variable_name

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
            # This is the first task. By definition, it has no dependencies
            task.mark_dependencies_found()
        else:
            for i in task.inputs:
                searched_task_index = task_index - 1
                dependency_found = False
                while searched_task_index >= 0 and not dependency_found:
                    searched_task = self.task_list[searched_task_index]

                    # Do not try to search delivery tasks, move up the task tree
                    if not searched_task.delivery_task:

                        # Looping through the outputs of the current task being searched
                        for output in searched_task.outputs:
                            # If the input/output filetypes match
                            if i.name == output.name:
                                task.add_formatted_input_string(
                                    Workflow.__format_import_string(i, output, searched_task)
                                )
                                dependency_found = True
                                break
                    searched_task_index -= 1

                # If this is reached, there are no more nodes that can be searched. The dependency could not be found
                if not dependency_found:
                    self.project_logger.log_error(
                        "E.wfg.InO.1",
                        "The " + i.name.upper() + " input within the " + task.alias.upper() +
                        " task could not find a fitting output in any of the upstream tasks"
                    )
                    sys.exit(1)
            task.mark_dependencies_found()

    def __process_task_list(self):
        """
        For each task, construct the input strings that the task definition will need

        side-effect of adding these input strings to the tasks internal list
        """
        for index, task in enumerate(self.task_list):
            self.__find_dependencies(task, index)

    def __construct_import_lines(self):
        """
        Construct the import lines needed by the WDL script

        :return: A collection of import lines as strings
        """
        import_lines = []
        for task in self.task_list:
            import_lines.append('import "' + task.import_location + '" as ' + task.alias.upper())
        return import_lines

    @staticmethod
    def __construct_task_lines(task: Task) -> List[str]:
        """
        Construct the task definition for a given task

        :param task: the input task
        :return: The task definition as a list of strings
        """
        lines = []
        alias = task.alias
        num_inputs = len(task.formatted_input_strings)

        first_line_fragment = Workflow.INDENT + "call " + alias.upper() + "." + task.task_name + " as " + alias.lower()

        if num_inputs > 0:
            lines.append(first_line_fragment + " {")
            lines.append(Workflow.INDENT * 2 + "input:")
            # Inputs with comma (all but the last statement)
            for i in task.formatted_input_strings[:-1]:
                lines.append(Workflow.INDENT * 3 + i + ",")
            # Last input (has no comma)
            lines.append(Workflow.INDENT * 3 + task.formatted_input_strings[-1])
            lines.append(Workflow.INDENT + "}")
        else:
            lines.append(first_line_fragment)
        return lines

    def __construct_output_lines(self):
        """
        Construct all of the lines that make up the WDL script

        :return: The lines of the WDL script
        """
        import_lines = self.__construct_import_lines()

        workflow_line = "workflow " + self.workflow_name + " {"

        task_lines = []
        for i in self.task_list:
            for lines in self.__construct_task_lines(i):
                task_lines.append(lines)
            task_lines + self.__construct_task_lines(i)
            task_lines.append("")

        final_line = "}"

        return import_lines + ["", workflow_line] + task_lines + [final_line]

    def create_and_write_wdl_file(self, output_file):
        """
        Create the WDL script and write it to a file

        :param output_file: the name of the output file
        """
        self.__process_task_list()
        lines = self.__construct_output_lines()

        with open(output_file, "w") as F:
            for line in lines:
                F.write(line + "\n")
