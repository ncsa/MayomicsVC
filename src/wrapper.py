#!/usr/bin/env python3

import sys
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('program', type=str, help='program to execute')
parser.add_argument('exit-on-failure', type=bool, help='True or False')

args = parser.parse_args()

def main(args):
    program = args.program
    exit = args.exit_on_failure

    result = subprocess.run(program, shell=True)

    if result.returncode != 0:
        if exit:
            # < Insert write to log file here >
            sys.exit(result.returncode)
        else:
            # < Figure out how to tell the next program not to run >
            # This is an "artificial zero", created just to keep cromwell from killing the workflow
            sys.exit(0)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
