# Purpose 

To purpose of the Pconvert the configuration files to a JSON format, we use a Java program called WOM tool. The reason being because JSON is the format used by Cromwell which interprets the WDL Programs.
However, when the JSON file is complete, it mentions the types of the variables and not the actual variables present in the JSON file. Hence, to convert the JSON file so that it mentions the actual variables, we use a python program named config_parser that takes the JSON format of all the variables with each of their types and the configuration files as an input and outputs a JSON file with variables filled in.

