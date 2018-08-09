#!/usr/bin/env python

###
#  All names in this file are lower case, because all checks of key names are made lower case with the lower() method
###

# This is the list of keys that, if turned off, must be set to the string "null" ('Key="null") in the flat config file
NULLABLE_KEYS = ("inputread2")

# This is the list of keys that are allowed to be blank (i.e. set to either 'Key=' or 'Key=""') in the flat config file
OPTIONAL_KEYS = ("debugmode")

