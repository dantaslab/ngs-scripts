#!/usr/bin/env python

"""
File Name    : utilities_MB.py
Author       : Manish Boolchandani, manish@wustl.edu
Created On   : 2017-04-10
Last Modified: 2016-07-05

Description  :  This utility module contains function that are used by
                multiple programs

Dependencies:

CHANGE LOG:
TODO:

"""

# Python imports
import sys
import os
import time
import re
import argparse
import logging
import subprocess

# Logger configuration and setup
logger = logging.getLogger(__name__)


def file_exists_readable(file, raise_IOError=None):
    """
    Exit with error if file does not exist or is not readable
    Or raise an IOerror if selected
    """

    if not os.path.isfile(file):
        message="Can not find file "+ file
        logger.critical(message)
        if raise_IOError:
            print("CRITICAL ERROR: " + message)
            raise IOError
        else:
            sys.exit("CRITICAL ERROR: " + message)

    if not os.access(file, os.R_OK):
        message="Not able to read file " + file
        logger.critical(message)
        if raise_IOError:
            print("CRITICAL ERROR: " + message)
            raise IOError
        else:
            sys.exit("CRITICAL ERROR: " + message)


def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """

    paths = os.environ["PATH"].split(os.pathsep)

    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return fullexe
    return


def check_outfiles(outfiles):
    """
    If outfiles already_exist, then remove or bypass
    """
    bypass=[]
    for file in outfiles:
        if os.path.isfile(file):
            if config.resume and os.path.getsize(file) > 0:
                bypass.append(True)
            else:
                bypass.append(False)
        else:
            bypass.append(False)

    if False in bypass or not bypass:
        # remove any existing files
        for file in outfiles:
            remove_file(file)
        return False
    else:
        return True

"""
# Parse step number argument and returns the list of steps to run in the pipeline
#
"""
def process_StepCount(step_num, MAXSTEP):

        procStep = []

        if step_num == "0":
                for i in range(1, MAXSTEP+1):
                        procStep.append(i)
        else:
                steps = step_num.split(':')

                if len(steps) == 2:
                        start = int(steps[0])
                        end = int(steps[1])

                        if end > MAXSTEP: end = MAXSTEP #       logger.warning("Max Step Value is {}".format(MAXSTEP))
                        if start <= 0: start = 1
                        if end > start:
                                for step in range(start, end + 1):
                                        procStep.append(step)
                        else:
                                logger.exception("Error: Step information not in correct order. Enter analyze_metagenome -h to check the usage")
                                sys.exit(1)
                elif len(steps) == 1 and 0 < int(steps[0]) <= MAXSTEP:
                        procStep.append(steps[0])

                else:
                        logger.exception("Error: Step information not in correct order. Enter analyze_metagenome -h to check the usage")
                        sys.exit(1)

        logger.info("Steps: {}".format(str(procStep)))

        return procStep

"""
# This function will submit jobs to cluster
#
"""
def run_sbatch_script(script_file, split_cnt, workdir, sbatch_params, temp_dir):

    logger.info("Submitting jobs to cluster")
    job_status = "FAIL"
    #sbatch_script = '/scratch/gdlab/manish/scripts/ngs/submit_job_array_slurm.sh'
    #sbatch_script = '/opt/apps/labs/gdlab/software/ngs-scripts/1.0/submit_job_array_slurm.sh'
    sbatch_script = 'submit_job_array_slurm.sh'
    #logger.debug("Inside function:" + script_file)

    job_status = "FAIL"
    jobId = None

    if not sbatch_params:
        sbatch_params = '--mem=16000M'

    if not split_cnt:
        split_cnt = 1

    cmd = "{} {} {} | sbatch {}".format(sbatch_script, script_file, str(split_cnt), str(sbatch_params));
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        logger.error("Error submitting script to cluster = {}".format(cmd))
        logger.error("Error Message ={}".format(str(e)))
        sys.exit(1)

    matchObj = re.search('^Submitted batch job (\d+)', output)

    if matchObj:
        jobId = matchObj.group(1)
        logger.debug("jobId = " + jobId)
    else:
        print "No match!"
        logger.exception("Error occured: No match!!. Could not get the job Id")
        sys.exit(1)

    if jobId:
        cmd2 = "sbatch --dependency=afterany:" + str(jobId) + " check_job_completion_MB.py \
                --job {} --workdir {} --tempdir {}".format(str(jobId), str(workdir), str(temp_dir))
	logger.debug("cmd2 = {}".format(str(cmd2)))

        try:
            filename = os.path.join(str(workdir), str(temp_dir), "check_" + str(jobId))
            logger.debug("filename = {}".format(str(filename)))

            cmd2_stat = subprocess.check_output(cmd2, shell=True, stderr=subprocess.STDOUT)
	    matchObj2 = re.search('^Submitted batch job (\d+)', cmd2_stat)

            if matchObj2:
	        temp_jobId = matchObj2.group(1)
        	logger.debug("temp jobId = " + temp_jobId)
    	    else:
        	print "No match!"
        	logger.exception("Error occured: No match!!. Could not get the job Id")
        	sys.exit(1)

            while not os.path.exists(filename):
                #logger.debug("Waiting for a job to finish!")
                time.sleep(5)

            if os.path.exists(filename):
                job_status = "SUCCESS"
		cmd3 = "mv slurm-{}* {}/".format(str(jobId), os.path.join(workdir, temp_dir))
		cmd3_stat = subprocess.check_output(cmd3, shell=True, stderr=subprocess.STDOUT)
		cmd4 = "mv slurm-{}* {}/".format(str(temp_jobId), os.path.join(workdir, temp_dir))
		cmd4_stat = subprocess.check_output(cmd4, shell=True, stderr=subprocess.STDOUT)

            else:
                logger.info("job failed!")
                sys.exit(1)

        except subprocess.CalledProcessError as e:
            logger.error("Error submitting script to cluster = {}".format(cmd2))
            logger.error("Error Message = {}".format(str(e)))
            sys.exit(1)

            #status = wait_for_completion(jobId)
    else:
        logger.exception("Invalid job number!")
        sys.exit(1)

    return job_status;
