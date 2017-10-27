1 Objective
============

Recreate GenomeGPS in Cromwell/WDL instead of Bash/Perl

2 Workflow architecture
=======================

2.1 Basic coding principles
---------------------------

* Nothing should be hard-coded - paths to executables, software parameters, GATK bundle files should be defined in the runfile
* Comments in code
* Documentation in readme files



2.2 Workflow best practices
---------------------------

* Must be modular to be maintainable
* Must be robust against hardware/software/data failure
    * user option on whether to fail or continue the whole workflow when something goes wrong with one of the sample
    * produce logs on failure; capture exit codes; write to FAIL log file; email analyst 
    * check everything before workflow actually runs: 
       * check that all the sample files exist and have nonzero size 
       * check that all executables exist
       * for each workflow module at runtime, check that output was actualy produced and has nonzero size
       * perform QC on each output file, write results into log, give user option to continue anyway if QC is failed
* Unit testing
    * Must have tests designed to succeed and designed to fail to ensure that failure mechanisms work




2.3 Need workflow diagram here
------------------------------

* Reconstruct from GenomeGPS file structure, highlight which parts we are doing in what order.
* Mayo and UIUC will have division of labor among the modules - highlight those too, in different color



3 Dependencies
==============

Grab from GenomeGPS documentation, highlight which parts we are doing in what order.

4 To-Dos
========

* As of Oct 24, 2017: Ram will work only on the bwa-mem module, to implement fully the template that could be used for other modules
    * multiple samples in parallel
    * all paths and parameters coded in runfile
    * loggery, user notification
    * error capture
    * checking inputs and executables
    * checking outputs
    * QC on input FASTQ
    * QC on outputs
```
  |
  |-call-Samtools
  |  |-shard-2
  |  |  |-inputs
  |  |  |  |-projects
  |  |  |  |  |-mgc
  |  |  |  |  |  |-Project_1
  |  |  |  |  |  |  |-ram
  |  |  |  |  |  |  |  |-CromwellWDL
  |  |  |  |  |  |  |  |  |-MultiSampleMultiStepVC
  |  |  |  |  |  |  |  |  |  |-cromwell-executions
  |  |  |  |  |  |  |  |  |  |  |-BWA_Mem_Run
  |  |  |  |  |  |  |  |  |  |  |  |-97cbc5de-ff36-4912-aa04-395b08702c85
  |  |  |  |  |  |  |  |  |  |  |  |  |-call-BWA_Mem
  |  |  |  |  |  |  |  |  |  |  |  |  |  |-shard-2
  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |-execution
  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |-glob-8fbfe7a84f921347caf0a4408578e296
  |  |  |-execution
  |  |  |  |-tmp.Jh4vIo
  |  |  |  |-glob-e04fb83eeaea0b4b8fdde6282d30506c
  |  |-shard-0
  |  |  |-inputs
  |  |  |  |-projects
  |  |  |  |  |-mgc
  |  |  |  |  |  |-Project_1
  |  |  |  |  |  |  |-ram
  |  |  |  |  |  |  |  |-CromwellWDL
  |  |  |  |  |  |  |  |  |-MultiSampleMultiStepVC
  |  |  |  |  |  |  |  |  |  |-cromwell-executions
  |  |  |  |  |  |  |  |  |  |  |-BWA_Mem_Run
  |  |  |  |  |  |  |  |  |  |  |  |-97cbc5de-ff36-4912-aa04-395b08702c85
  |  |  |  |  |  |  |  |  |  |  |  |  |-call-BWA_Mem
  |  |  |  |  |  |  |  |  |  |  |  |  |  |-shard-0
  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |-execution
  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |-glob-e2984cba736e7acdf7c0536378232fca
  |  |  |-execution
  |  |  |  |-tmp.JOg8Ce
  |  |  |  |-glob-e04fb83eeaea0b4b8fdde6282d30506c
  |  |-shard-1
  |  |  |-inputs
  |  |  |  |-projects
  |  |  |  |  |-mgc
  |  |  |  |  |  |-Project_1
  |  |  |  |  |  |  |-ram
  |  |  |  |  |  |  |  |-CromwellWDL
  |  |  |  |  |  |  |  |  |-MultiSampleMultiStepVC
  |  |  |  |  |  |  |  |  |  |-cromwell-executions
  |  |  |  |  |  |  |  |  |  |  |-BWA_Mem_Run
  |  |  |  |  |  |  |  |  |  |  |  |-97cbc5de-ff36-4912-aa04-395b08702c85
  |  |  |  |  |  |  |  |  |  |  |  |  |-call-BWA_Mem
  |  |  |  |  |  |  |  |  |  |  |  |  |  |-shard-1
  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |-execution
  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |-glob-1b331778c20fad9525cbbad0d5f04486
  |  |  |-execution
  |  |  |  |-glob-e04fb83eeaea0b4b8fdde6282d30506c
  |  |  |  |-tmp.yCh857
  |-call-BWA_Mem
  |  |-shard-2
  |  |  |-inputs
  |  |  |  |-projects
  |  |  |  |  |-bioinformatics
  |  |  |  |  |  |-DataPacks
  |  |  |  |  |  |  |-human
  |  |  |  |  |  |  |  |-EllensSyntheticData
  |  |  |  |  |  |  |  |  |-WGS_chr1_5X_E0.005
  |  |  |  |  |  |  |  |  |  |-Chunks
  |  |  |  |  |  |-TestData
  |  |  |  |  |  |  |-HG19
  |  |  |  |  |  |  |  |-gatk_bundle_2.8
  |  |  |  |  |  |  |  |  |-2.8
  |  |  |  |  |  |  |  |  |  |-hg19
  |  |  |-execution
  |  |  |  |-glob-8fbfe7a84f921347caf0a4408578e296
  |  |  |  |-tmp.OWyZa0
  |  |-shard-0
  |  |  |-inputs
  |  |  |  |-projects
  |  |  |  |  |-bioinformatics
  |  |  |  |  |  |-DataPacks
  |  |  |  |  |  |  |-human
  |  |  |  |  |  |  |  |-EllensSyntheticData
  |  |  |  |  |  |  |  |  |-WGS_chr1_5X_E0.005
  |  |  |  |  |  |  |  |  |  |-Chunks
  |  |  |  |  |  |-TestData
  |  |  |  |  |  |  |-HG19
  |  |  |  |  |  |  |  |-gatk_bundle_2.8
  |  |  |  |  |  |  |  |  |-2.8
  |  |  |  |  |  |  |  |  |  |-hg19
  |  |  |-execution
  |  |  |  |-glob-e2984cba736e7acdf7c0536378232fca
  |  |  |  |-tmp.rrJiKx
  |  |-shard-1
  |  |  |-inputs
  |  |  |  |-projects
  |  |  |  |  |-bioinformatics
  |  |  |  |  |  |-DataPacks
  |  |  |  |  |  |  |-human
  |  |  |  |  |  |  |  |-EllensSyntheticData
  |  |  |  |  |  |  |  |  |-WGS_chr1_5X_E0.005
  |  |  |  |  |  |  |  |  |  |-Chunks
  |  |  |  |  |  |-TestData
  |  |  |  |  |  |  |-HG19
  |  |  |  |  |  |  |  |-gatk_bundle_2.8
  |  |  |  |  |  |  |  |  |-2.8
  |  |  |  |  |  |  |  |  |  |-hg19
  |  |  |-execution
  |  |  |  |-tmp.kbVcdk
  |  |  |  |-glob-1b331778c20fad9525cbbad0d5f04486
```
