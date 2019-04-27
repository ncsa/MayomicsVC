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
 
