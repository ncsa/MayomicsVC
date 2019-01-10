# Testing the Shell Scripts

## Assumed folder structure

These tests as written assume you test from your cloned repo from this folder. The python code shell_script_testing.py is intended to be run from this folder. Assuming the repo is cloned in /home/username/MyTestingSpree, you'll need the following folders in that directory:
* MayomicsVC/ - the actua repo
  * DON'T PUT ANYTHING IN HERE ITHAT IS NOT PART OF THE ACTUAL CODE BASE
* Config/ You'll need an EnvProfile.file to run the tests as is. This file is simply a file with the Sentieon license
* Inputs/ This folder will contain any inputs files: fastq, bams, vcfs, as well as the TruSeqAdaptors.fast adaptor file. 
* Delivery/ This is where the deliver_XXX scripts will deliver the final files from the workflow
* Reference/ This will contain reference files including the human genome assembly reference and variant references
* Jsons/ This will contain the json file that the wdl workflow manager needs in order to run, which is delivered by the deliver_XXX scripts so that the workflow can be run again. For testing purposes, it need only be non-empty

## Executing the code

To execute the tests as written, you'll need Python3 installed, as well as python's cutadapt tool. For many steps, you'll also need the Sentieon tools installed and a license to run them. Check the tests for expected input format for each script. To run the test, from a Python environment, you'll run the code with the intended script to test as an argument.

```shell
>>> python shell_script_testing.py trim_sequences.sh
```
