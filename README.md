1 Objective
============

Recreate GenomeGPS in Cromwell/WDL instead of Bash/Perl

2 Workflow architecture
=======================

1.2 Basic coding principles
---------------------------

1. Nothing should be hard-coded - paths to executables, software parameters, GATK bundle files should be defined in the runfile

2. Comments in code

3. Documentation in readme files



1.2 Workflow best practices
---------------------------

1. Must be modular to be maintainable

2. Must be robust against hardware/software/data failure

    2.1 user option on whether to fail or continue the whole workflow when something goes wrong with one of the sample

    2.2 produce logs on failure; capture exit codes; write to FAIL log file; email analyst 

    2.3 check everything before workflow actually runs: 

       * check that all the sample files exist and have nonzero size 

       * check that all executables exist

       * for each workflow module at runtime, check that output was actualy produced and has nonzero size

       * perform QC on each output file, write results into log, give user option to continue anyway if QC is failed

    2.4 Must have tests designed to succeed and designed to fail to ensure that failure mechanisms work


3 Dependencies
==============
