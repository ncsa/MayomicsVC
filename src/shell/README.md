# Purpose
 
The purpose of this readme is to understand the peculiarities of the bash scripts. The shell scripts are structured in a way that can be distinguished into two sections, standard and variable. The standard section consists of the lines that are similar across all bash scripts except for file name change, whereas, the variable section is individual to each bash script.
 
# Bash Scripting Structure
 
<img src="https://user-images.githubusercontent.com/43070131/56844379-0708f280-6875-11e9-9efa-7b03910e5917.PNG" width="800">

# Description

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
