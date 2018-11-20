# Testing procedures for MayomicsVC

# Assumed folder structure

Wherever you test, i.e. /home/username/MyTestingSpree/, the following folder structure helps to keep it all tidy:
* Config/ - will have your environmental profile files and the config files
* Jsons/ - where the JSON files will go, after WOMTOOL creates them and parser parses them
* Inputs/ - where the reference and input fastq files reside - but it's also ok to point config files to another, central location, accessible to everyone
* Delivery/ - where the output files from the workflow will be placed by the Delivery blocks
* MayomicsVC/ - the actual repo
  * DON'T PUT ANYTHING IN HERE THAT IS NOT PART OF THE ACTUAL WDL/BASH/PYTHON_PARSER CODE BASE
* cromwell-executions/ - folder created by Cromwell when running the workflow; this is where debugging happens
* cromwell-workflow-logs/ 

# The vision

Ideally we will grow a battery of feature tests for Bash and separately a battery of integration tests for WDL.

For example, right now there is only RunWdlIntegrationTest.sh in this space. It will run a single integration test on one workflow specified as input to that shell script. We may want to run multiple of these at the push of a button. Thus we may script a series of these, and start them from within a single script. This could be automated via Jenkins.

Similarly, we are writing Bash integration tests in Bash. These will become part oif this subfolder and be invoked from the same single script.
