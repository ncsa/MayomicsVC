# Purpose 

To convert the configuration files to a JSON format, we use a Java program called WOM tool. The reason being because JSON is the format used by Cromwell which interprets the WDL Programs.
However, when the JSON file is complete, it mentions the types of the variables and not the actual variables values present in the file. The user fills the JSON values manually. Hence, to automatically fill the JSON file with the variable values, we use a python program named config_parser that takes the JSON format of all the variables with each of their types and the configuration files as an input and outputs a JSON file with variables filled in.

