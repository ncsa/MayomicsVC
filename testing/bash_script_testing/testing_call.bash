#!/bin/bash

export TEST_PATH="$BATS_TEST_DIRNAME"
export MAYOMICS_PATH="$(cd "$TEST_PATH/../.." && pwd)"
export SCRIPT_PATH="$MAYOMICS_PATH/src/bash_scripts"
export SCRIPT_TO_TEST=$1
