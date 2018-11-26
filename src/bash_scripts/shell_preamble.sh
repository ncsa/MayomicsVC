#!/bin/bash

set -euxo pipefail

function sigusrhandler1()
{
    echo "SIGUSR1 caught by shell script" 1>&2
    echo 30 > ./rc
    sync
}

function sigusrhandler2()
{
    echo "SIGUSR2 caught by shell script" 1>&2
    echo 31 > ./rc
    sync
}

trap sigusrhandler1 SIGUSR1
trap sigusrhandler2 SIGUSR2


export LD_LIBRARY_PATH=""