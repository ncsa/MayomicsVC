# Background

We need a way to validate that all of the Cromwell/WDL workflow configuration variables are valid:

For example, if a variable is supposed to point to an executable, we should programmatically
  confirm that the file exists and that one has permission to execute it.

## Problem

The JSON config template file that Cromwell/WDL generates has types, but they are not specific enough
  i.e. executables are simply listed as "Strings", only inputs are listed as "Files", etc.

We need more specific typing information in order to know programmatically how to validate that a given key is pointing
  to an valid value

## Solution

We construct a json file that contains the typing information of all of the variables present in the workflow
  configuration files.

The 'key_types.json' file contains only one kind of entry: key-value pairs.
  Where the key represents a variable present in the Cromwell/WDL config file
  and the value represents the keys type

### Valid Types

#### Arrays of Types
    
    For any arrays in the type file, write square brackets '[]' around the type.
    Even if the array is nested, it only requires single brackets.
    All lists entries must be of the same type
    
    Examples:
        For an array of files, write "[Files]"
        For an array of array of Integers, write "[Int]"  

#### ExecFile

For ExecFile types, we confirm that
    
    1. The value is a file that exists
    2. The caller has permission to execute the file

#### File

For File types we confirm that the value is a file that exists

#### Boolean

For Boolean types we confirm that the value EXACTLY matches either 'true' or 'false'

This was made case sensitive due to how the BASH scripts were written

#### String

The String type is to generic to check during the pre-workflow quality control.
They are just passed without any checks.

#### Integer

Integer types are validated by seeing if python can convert them to an int. 
  Values like 7 or 40000 will pass, but values like 3.14159 or "Eight" will fail

#### Decimal

Decimal types are validated by seeing if python can convert them to a floats
    
#### Directory

For this basic Directory type, we confirm that the directory exists
  and do not check any permissions.
  
If one wants read/write/executable permissions checks, use a more
  specific directory type below.
    
#### OutputDirectory

For types listed as OutputDirectory, we confirm that

    1. The directory does not already exist

       
This type should be used for directories that the workflow will write into.

#### ReadOnlyDirectory

For types listed as ReadOnlyDirectory, we confirm that

    1. The directory exists
    2. It has readable and executable permissions
       for the current user
       
Readable means that the user can read the files inside of the directory, and
executable means that the user can move into the directory.
       
This type should only be used for directories that only need to be read.

#### DebugMode

This is a special case. The variable 'DebugMode' has only one acceptable value: '-d'.
  This value is passed directly to the bash scripts, and passing any other value (other
  than an empty value) will kill the script.

# IMPORTANT NOTE!

 If any new variables are added to the Cromwell/WDL code, 
 their types must be added to the types list, located at:
 
 ```src/config/key_types.json```
