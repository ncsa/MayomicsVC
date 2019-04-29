# Purpose
 
The purpose of this readme is to understand the peculiarities of the bash scripts. The bash scripts are structured in a way that can be distinguished into two sections, standard and variable. The standard section consists of the lines that are similar across all bash scripts except for file name change, whereas, the variable section is individual to each bash script. The script used to explain the bash structure is <a href="https://github.com/ncsa/MayomicsVC/blob/master/src/shell/trim_sequences.sh"> Trim_Sequences </a>
 
# Bash Scripting Structure
 
<img src="https://user-images.githubusercontent.com/43070131/56844379-0708f280-6875-11e9-9efa-7b03910e5917.PNG" width="800">

# Description

## Standard

<details>
  <summary>
  Manifest, Usage Docs and Set Checks
  </summary> 
 
 ```bash scripting
 #!/bin/bash
 ```
The #!/bin/bash command is mandatory and it is best practices for script writing in bash.
It is a standalone executable and when invoked, ./filename.sh needs to be typed and it should run.

 ```bash scripting
read -r -d '' MANIFEST << MANIFEST
 ```
MANIFEST is a built-in bash command. This command is essential for debugging because if the workflow was executed and variants were called, users should be able to track the history. 
 
 ```bash scripting
 *****************************************************************************
`readlink -m $0`
called by: `whoami` on `date`
command line input: ${@}
*****************************************************************************
 ```
This command is essential for debugging because it stores the details about who ran the workflow and when.

```bash scripting
set -o errexit
```
This statement requires that bash scripts exit whenever there is an error.
Usually, bash scripts continue to run even if one command failed.
-o errexit will prevent that from happening and will quit if any error occurs at all.

```bash scripting
set -o pipefail
```
-o pipefail sets the exit code in the pipeline to the right most command to exit with a non-zero status.

```bash scripting
set -o nounset
```
The -o nounset is used to error out any variable that is not defined in the script.
This requires that all variables are defined and it should not have any loose variables that are no longer needed.

```bash scripting
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
 
 ```bash scripting
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
 
 ```bash scripting
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

## Variable
<details>
<summary>
Precheck for Input and Output
</summary> 
 
```bash scripting
source ${SHARED_FUNCTIONS}
 
checkVar "${SAMPLE+x}" "Missing sample name option: -s" $LINENO

ERRLOG=${SAMPLE}.trimming.${SGE_JOB_ID}.log
    ....
    ....
echo "${MANIFEST}" >> "${ERRLOG}"
  
checkVar "${ENV_PROFILE+x}" "Missing environmental profile option: -e" $LINENO
source ${ENV_PROFILE}
  
checkVar "${ADAPTERS+x}" "Missing adapters file option: -A" $LINENO
        ....
        ....
checkVar "${INPUT2+x}" "Missing read 2 option: -r. If running a single-end job, set -r null in command." $LINENO
  
checkVar "${IS_PAIRED_END+x}" "Missing paired-end option: -P" $LINENO
  
if [[ "${IS_PAIRED_END}" != true ]] && [[ "${IS_PAIRED_END}" != false ]]
then
        ....
        ....
fi
if [[ "${IS_PAIRED_END}" == true ]]
then
        ....
        ....
fi
if [[ "${IS_PAIRED_END}" == false ]]
        ....
        ....
fi
 
checkVar "${CUTADAPT+x}" "Missing CutAdapt software path option: -C" $LINENO
checkDir ${CUTADAPT} "Cutadapt directory ${CUTADAPT} is not a directory or does not exist." $LINENO
checkVar "${THR+x}" "Missing threads option: -t" $LINENO
```
Precheck calls functions from the shared functions.sh file to perform the following operations:

1. Checks if the sample name variable exists or not
2. Creates log for JOB_ID/script and tool
3. Sends Manifest to the Log
4. Sources the file with environmental profile variables
5. Check if input files, directories, and variables are non-zero

In the case of adapters, if the adapters string is defined to the full path of the file than the variable is set and we do not need to check that file as it would have been checked by the parser.So, the argument to the full path to the adapter file + x ( "${ADAPTERS+x} ") will be passed into the checkVar.It will check it as the first variable of the string and will not throw an error by setting the error code not equal to 1.However, if the full path to the adapter file is not defined than string + x is passed and bash will pass the empty string as the first variable. The exit code will be set to 1 and an error will be thrown.
</details>

<details>
<summary>
FileName Parsing
</summary> 

```bash scripting
## Parse filename without full path
OUT1=$(basename ${INPUT1})
if  [[ "${IS_PAIRED_END}" == false ]]  # If single-end, we do not need a second output trimmed read
then
               OUT2=null
else
               OUT2=$(basename ${INPUT2})
fi
```

The filename parsing section parses the filename without the full path and the reason why it does without the full path is because  cutadapt requires the output option, -o, and hence, file name parsing is necessary.
</details>

<details>
 <summary>
  CutAdapt Read Trimming
 </summary>
 
 ``` Bash Scripting
 ## Record start time
logInfo "[CUTADAPT] START."
## Cutadapt command, run for each fastq and each adapter sequence in the adapter FASTA file.
## Allocates half of the available threads to each process.
if [[ "${IS_PAIRED_END}" == false ]]  # if single-end reads file
then
               # Trim single-end reads
               TRAP_LINE=$(($LINENO + 1))
               trap 'logError " $0 stopped at line ${TRAP_LINE}. Cutadapt Read 1 failure. " ' INT TERM EXIT
               ....
               ....
else
               TRAP_LINE=$(($LINENO + 1))
               trap 'logError " $0 stopped at line ${TRAP_LINE}. Cutadapt Read 1 and 2 failure. " ' INT TERM EXIT
               ....
               ....
 
        checkExitcode ${EXITCODE} $LINENO
               logInfo "[CUTADAPT] Trimmed adapters in ${ADAPTERS} from input sequences. CUTADAPT log: ${SAMPLE}.cutadapt.log"
```

The Cutadapt command runs for each fastq and each adapter sequence in the adapter fasta file. It allocates half of the available threads to each process. If the file is single-end reads than it trims the single end reads file. If the file is in paired-end mode, than trimming reads with cutadapt occurs in paired-end mode, where -a and -A specify forward and reverse adapters, respectively. -p specifies output for read2.

The reason why traps are used in the cutadapt read trimming is mentioned below:

The command, set -o error exit is mentioned because the script should be forced to exit whenever there is an issue. However, the error statement would not record the error to the right log and would not report the right line number. The exit error statement would force the bash to quit at the cutadapt line. So, the trap is specified in a nested way where the outer trap contains the statement that the error log would print whenever there is an interruption, termination or an exit in the workflow. The important components of traps are:

1. Trap_line : Variable that refers to the line number where the error occurred
2. Log error : Function needs to be inside the quotes for the trap to act on it
3. Dollar zero ($0) : Refers to the name of the bash script followed by input parameters
4. checkExitcode : checkExitcode checks whether the exit code is zero or not. It needs only two inputs, So it only needs 2 inputs, the exit code and the line number because the function needs to print the exit code with the message that is meaningful into the log and it needs to state the line number at which the exit code is non-zero.
</details>               


<details>
 <summary>
Post Processing
 </summary>
 
 ```bash scripting
 ## Check for file creation
checkFile ${OUT1} "Output trimmed read 1 file ${OUT1} is empty." ${LINENO}
if [[ "${IS_PAIRED_END}" == true ]]
then
        checkFile ${OUT2} "Output trimmed read 2 file ${OUT2} is empty." ${LINENO}
fi
```
The post processing section checks for file creation and if the read files are empty or not
               
               
