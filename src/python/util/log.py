#!/usr/bin/env python
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 6:
    print("The script (" + sys.argv[0] + ") requires Python version 3.6 or higher")
    sys.exit(1)

import logging

'''
Instead of using Python's logger module directly, this class is defined to abstract the details away from the client
  code

Mayo has standardized on the following log levels: INFO, DEBUG, WARN, and ERROR

Log message format:
  TIMESTAMP(YYYY-MM-DD'T'HH:mm:ssz) SEVERITY SCRIPTNAME JOBID/NA ERRORCODE MESSAGE
  Example: 2016-03-25T13:01:22-500 ERROR myScript.sh 929122 E32912 input bam file does not have sufficient reads

The JOBID/NA and ERRORCODE terms cannot be easily added to the default logging format, and must be included in
  the logged message directly. To make this easy, this class defines wrapper methods around calls to the logger methods
  that automatically format the JOBID and ERRORCODE terms
'''
class ProjectLogger:
    def __init__(self, job_id, name, log_level=logging.INFO):
        logging.basicConfig(level=log_level,
                            format="[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S%z"
                            )
        self.logger = logging.getLogger(name)
        self.job_id = str(job_id)
        self.warnings_issued = 0
        self.errors_issued = 0

    def log_info(self, message):
        self.logger.info("[" + self.job_id + "] [-] " + message)

    def log_debug(self, message):
        self.logger.debug("[" + self.job_id + "] [-] " + message)

    def log_warning(self, message):
        self.logger.warning("[" + self.job_id + "] [-] " + message)
        self.warnings_issued += 1

    def log_error(self, error_code, message):
        self.logger.error("[" + self.job_id + "] [" + error_code + "] " + message)
        self.errors_issued += 1
