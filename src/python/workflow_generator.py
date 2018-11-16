#!/usr/bin/env python3

from generators.workflows import Workflow
import sys
def main(args):

    workflow = Workflow("TestWorkflow", "generators/generator_test_resources/test_workflow.txt", False, job_id="NA")

    workflow.process_task_list()

    for i in workflow.construct_output_lines():
        print(i)# + "\n")




if __name__ == "__main__":
    main(sys.argv[1:])