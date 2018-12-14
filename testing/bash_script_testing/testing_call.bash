#!/bin/bash

export TEST_PATH="$BATS_TEST_DIRNAME"
export Mayomics_path="$(cd "$TEST_PATH/../.." && pwd)"
export SCRIPT_PATH="$Mayomics_path/src/bash_scripts"
export SCRIPT_TO_TEST=$1
export Input_File_Path="$BATS_TEST_DIRNAME/../../.."
