#!/usr/bin/env python3

from typing import List
import collections
from generators.file_types import *
from generators.tasks import *


sorted_task_list: List[Task] = sorted(task_dict.values(), key=lambda task: task.rank)
sorted_task_names: List[str] = [i.alias.upper() for i in sorted_task_list]


def __get_start_list():
    """
    The full list of possible inputs at the beginning of the workflow building process (leaves off the last task)
    """
    return [i.alias.upper() for i in sorted_task_list[:-1]]


def __create_prompt(prompt_text, task_choices: List[str]) -> str:
    """
    Creates the interactive prompt string, and returns the input passed in from the end user
    Will fail if the user passes in invalid input over 5 times
    """
    fail_count = 5
    full_text = prompt_text + "[" + " | ".join(task_choices) + "]\n>>>"

    while fail_count > 0:
        input_string = input(full_text).upper()
        if input_string in task_choices:
            return input_string
        else:
            print("Error: " + input_string + " is not a valid input")
            fail_count -= 1
    print("Input invalid input passed in too many times. Exiting")


def __find_valid_next_choices(task_index) -> List[str]:
    """
    Given the index of a task, list all of the downstream tasks whose inputs match the current tasks outputs up to (and
     possibly including) a required task

    :param task_index: the index of the task to find the downstream choices for
    :return: list of valid next task choices
    """
    valid_choices = []
    task_outputs: List[FileType] = sorted_task_list[task_index].outputs

    curr_index = task_index + 1
    while curr_index < len(sorted_task_names):
        curr_task = sorted_task_list[curr_index]
        curr_task_name = sorted_task_names[curr_index]
        curr_index += 1

        # If the inputs and outputs between tasks match
        if collections.Counter(task_outputs) == collections.Counter(curr_task.inputs):
            valid_choices.append(curr_task_name)
        if curr_task.required:
            # Return because a required task was found
            return valid_choices
    # Return because the last task has been reached
    return valid_choices


def generate_task_list_interactively():
    workflow_tasks = []

    # Choose start point
    starting_choice = __create_prompt("Select the initial task: ", __get_start_list())
    starting_choice_index = sorted_task_names.index(starting_choice)

    workflow_tasks.append(sorted_task_list[starting_choice_index])

    current_task_index = starting_choice_index

    while current_task_index < (len(sorted_task_list) - 1):
        choices = __find_valid_next_choices(current_task_index)
        choices.append("END")

        input_task_name = __create_prompt("Select the next task: ", choices)

        if input_task_name == "END":
            break
        else:
            current_task_index = sorted_task_names.index(input_task_name)
            # Add the input to the workflow list
            workflow_tasks.append(sorted_task_list[current_task_index])

    return workflow_tasks
