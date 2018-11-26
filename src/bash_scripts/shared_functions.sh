#!/bin/bash

#########################################################
#
#  Logging functions used in MayomicsVC bash scripts
#
#########################################################

#-------------------------------------------------------------------------------------------------------------------------------

# Get date and time information
function getDate()
{
    echo "$(date +%Y-%m-%d'T'%H:%M:%S%z)"
}

# This is "private" function called by the other logging functions, don't call it directly,
# use logError, logWarn, etc.
function _logMsg () {
    echo -e "${1}"

    if [[ -n ${ERRLOG-x} ]]; then
        echo -e "${1}" | sed -r 's/\\n//'  >> "${ERRLOG}"
    fi
}

function logError()
{
    local LEVEL="ERROR"
    local CODE="-1"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    >&2 _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"

    if [[ -z ${EXITCODE+x} ]]; then
        EXITCODE=1
    fi

    exit ${EXITCODE};
}

function logWarn()
{
    local LEVEL="WARN"
    local CODE="0"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"
}

function logInfo()
{
    local LEVEL="INFO"
    local CODE="0"

    if [[ ! -z ${2+x} ]]; then
        CODE="${2}"
    fi

    _logMsg "[$(getDate)] ["${LEVEL}"] [${SCRIPT_NAME}] [${SGE_JOB_ID-NOJOB}] [${SGE_TASK_ID-NOTASK}] [${CODE}] \t${1}"
}

function checkArg()
{
    if [[ "${OPTARG}" == -* ]]; then
        echo -e "\nError with option -${OPT} in command. Option passed incorrectly or without argument.\n"
        echo -e "\n${DOCS}\n"
        exit 1;
    fi
}

#-------------------------------------------------------------------------------------------------------------------------------





##################################################################
#
# Functions to check for set variables and files 
#
##################################################################

#---------------------------------------------------------------------------------------------------------------------------------
#we need to figure out a good way to produce error messages, I think we should pass in a REASON string to argument $2, that's how I setting up the script for now.
#check for a set variale
function checkVar()
{
	if [[ ! -z ${1+x} ]]
	then
		EXITCODE=1
		logError "$0 stopped at line ${LINENO}. \nREASON=$2"
	fi

}


#check for the existence of a directory
#pass in the reason as the second parameter
function checkDir()
{
	if [[ ! -d $1 ]]
	then
		EXITCODE=1
		logError "$0 stopped at line $LINENO. \nREASON=$2"
	fi

}


#check for the existence of a file
#pass in the reason as the second parameter
function checkFile()
{
	if [[ ! -s $1 ]]
	then
        	EXITCODE=1
        	logError "$0 stopped at line $LINENO. \nREASON=$2"
	fi
}



