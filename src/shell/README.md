# Purpose
 
The purpose of this readme is to understand the peculiarities of the bash scripts. The shell scripts are structured in a way that can be distinguished into two sections, standard and variable. The standard section consists of the lines that are similar across all bash scripts except for file name change, whereas, the variable section is individual to each bash script.
 
# Bash Scripting Structure
 
<img src="https://user-images.githubusercontent.com/43070131/56844379-0708f280-6875-11e9-9efa-7b03910e5917.PNG" width="800">

# Description

<details>
<summary>
Standard
</summary>

<details>
  <summary>
   Manifest, Usage Docs and Set Checks
  </summary> 
 
 ```
 #!/bin/bash
 ```
The #!/bin/bash command is mandatory and it is best practices for script writing in bash.
It is a standalone executable and when invoked, ./filename.sh needs to be typed and it should run.

 ```
read -r -d '' MANIFEST << MANIFEST
 ```
MANIFEST is a built-in bash command. This command is essential for debugging because if the workflow was executed and variants were called, users should be able to track the history. 
 
 ```
 *****************************************************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*****************************************************************************
 ```
This command is essential for debugging because it stores the details about who ran the workflow and when.

```
set -o errexit
```
This statement requires that bash scripts exit whenever there is an error.
Usually, bash scripts continue to run even if one command failed.
-o errexit will prevent that from happening and will quit if any error occurs at all.

```
set -o pipefail
```
-o pipefail sets the exit code in the pipeline to the right most command to exit with a non-zero status.

```
set -o nounset
```
The -o nounset is used to error out any variable that is not defined in the script.
This requires that all variables are defined and it should not have any loose variables that are no longer needed.

```
SCRIPT_NAME=trim_sequences.sh
SGE_JOB_ID=TBD   # placeholder until we parse job ID
SGE_TASK_ID=TBD  # placeholder until we parse task ID
```
This script records the actual Job ID’s of every bioinformatics tasks that runs on the cluster.

</details>

<details>
  <summary>
   Logging Functions
  </summary>
 
 ```
 function checkArg()
{
    if [[ "${OPTARG}" == -* ]]; then
        echo -e "\nError with option -${OPT} in command. Option passed incorrectly or without argument.\n"
        echo -e "\n${DOCS}\n"
        exit 1;
    fi
}
 ```
The function checkArg() checks whether the option argument was passed correctly or not.
If any option is passed incorrectly, than the script will display the error "Error with option -${OPT} in command. Option passed incorrectly or without argument.\n"

</details>

<details>
  <summary>
   Getopts Argument Parser
  </summary>
 
 ```
 while getopts ":hl:r:A:C:t:P:s:e:F:d" OPT
do
               case ${OPT} in
                              h )  echo -e "\n${DOCS}\n"
                                             exit 0
                                         ...
                                         ...
                                         ...
                           : )     echo -e "\nOption -${OPTARG} requires an argument.\n\n${DOCS}\n" exit 1
                                     
               esac
done
 ```
 Principles for the Getopts Argument Parser of the script:

1. A colon after the letter means that it is mandatory and if the colon is not present, it means that it is optional.

2. Hence, we prepend the list before the colon because if an invalid option is provided, than the " \? " will be called. The only reason to allow the getops fuction to land to the "/?" case is if we have invalid option to prepend the list by a colon.

3. Each colon is being read separately. The getopts loop is reading consecutively. The case command assigns each argument entered to a variable and checks to make sure that a valid argument was entered for the options that require one. For example, ‘-d’ is the debug command and should never receive an argument following it. If it did, this would throw an error.

4. The colon at the beginning of the list turns off bash’s built-in error reporting, allowing us to handle errors with our checkArg function and the functions that follow, and allowing us to handle no arguments and bad arguments in a more meaningful way. If you pass in an option that is not recognized, the case statement will reach “/?” and it will print an invalid option statement. The final colon case is to ensure that every required option received an argument.
 
</details>
</details>
