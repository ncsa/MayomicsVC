#!/usr/bin/env python3

from generators.workflows import Workflow
from generators.interactive import generate_task_list_interactively
from generators.tasks import *
import sys
import argparse


def parse_args(args):
    """
    By default, argparse treats all arguments that begin with '-' or '--' as optional in the help menu
      (preferring to have required arguments be positional).

    To get around this, we must define a required group to contain the required arguments
      This will cause the help menu to be displayed correctly
    """
    parser = argparse.ArgumentParser(description="Either pass file in with '-i' or use '--interactive'")

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("-o", required=True, metavar='<File.wdl>', help='The location of the output wdl file')
    required_group.add_argument("--name", required=True, metavar='<String>', help='The name of the workflow')
    # Truly optional argument
    parser.add_argument("-i", required=False, metavar='<File>', help="A workflow task list file")
    parser.add_argument('--interactive', action="store_true", help='Interactive mode',
                        default=False, required=False
                        )
    parser.add_argument('--jobID', type=str, metavar='<String>', help='The job ID', default='NA', required=False)
    # Debug mode is on when the flag is present and is false by default
    parser.add_argument("-d", action="store_true", help="Turns on debug mode", default=False, required=False)
    return parser.parse_args(args)


def main(args):
    parsed_args = parse_args(args)

    if parsed_args.interactive:

        task_list = generate_task_list_interactively()

        # Print a summary message to stdOut
        print("Attempting to build a workflow with the following inputs\n")
        for task in task_list:
            print(task.get_input_output_summary_string())

        generator = Workflow(workflow_name=parsed_args.name,
                             debug_mode=parsed_args.d,
                             job_id=parsed_args.jobID,
                             input_file=None,
                             task_list=task_list
                             )
    else:
        generator = Workflow(workflow_name=parsed_args.name,
                             debug_mode=parsed_args.d,
                             job_id=parsed_args.jobID,
                             input_file=parsed_args.i,
                             task_list=None
                             )
    generator.create_and_write_wdl_file(parsed_args.o)


if __name__ == "__main__":
    main(sys.argv[1:])
