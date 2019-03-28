#!/usr/bin/env python

"""
File Name    : preprocess_seqreads_MB.py
Author       : Manish Boolchandani, manish@wustl.edu
Created On   : 2017-04-04
Last Modified: 2017-04-11

Description  :  Performs preprocessing on raw sequencing reads and
        gives the clean reads after quality trimming and removing host contaminants.
        Below are the steps that pipeline runs:

        Step 1: Counts the number of reads for each sample
        Step 2: Run fastqc on each sample and combine results using multiqc
        Step 3: Run trimmomatic using used defined parameters. Modify config
        file 'preprocess_config' to change the parameter settings
        Step 4: Run deconseq to remove human reads.
        Step 5: Fix paired end reads
        Step 6: Run fastqc again on cleaned sequencing reads for each sample
        Step 7: Counts the number of clean reads for each sample

        Final step : Generate summary for each samples

Dependencies:   fastqc, multiqc, bbtools, trimmomatic, deconseq

Usage: preprocess_reads.py [-h] [--version] -i INPUT_FILE -d SEQ_DIR -o OUTPUT_DIR --step STEP_NUM
                           [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
CHANGE LOG:
    v0.9.0: 1) added 7th step to count clean reads
	    2) created tmp directory to save slurm output files
    	    3) added multiqc output with fastqc output
    v0.6.2: added an input option for config file
    v0.6.1: added a parameter (sbatch_params) to run_sbatch_script() function.
        Now parsing sbatch memmory parameters specific to each step.
TODO:
    1) create a histogram of read distribtion
    2) add a paramter to overwrite files or append when rerunning the pipeline
    3) Remove temporary files and directory
"""

import sys

import os
import operator
import time
import re

import logging
import argparse
import subprocess
import ConfigParser
import tempfile
from collections import defaultdict

# Other imports
import preprocess_config_MB as config
import utilities_MB as util

# Version
VERSION = "0.9.0"

# Logger configuration and setup
logger = logging.getLogger(__name__)

# Maximum steps in the pipeline
MAX_STEPS = 7

config_opts = None


def main(argv):

    logger.info("Started Preprocessing Pipeline v" + VERSION)

    # Parse arguments from command line
    args = parse_arguments(argv)
    # print(args)

    # Check the required files, software and dependencies
    check_requirements(args)

    # Read config file and write config settings to the log file
    global config_opts
    config_opts = read_configuration(args)

    # print(config_opts.trimmomatic)
    # config.log_settings()

    # Check Step Number
    if not args.step_num:
        sys.exit("CRITICAL ERROR: Step number is not defined")
    else:
        steps_to_run = util.process_StepCount(args.step_num, MAX_STEPS)

    # step_info = dict.fromkeys(steps_to_run, 0)
    # logger.info("step_info = {}".format(str(step_info)))

    output_maindir = args.output_dir
    input_seqdir = args.seq_dir
    input_file = args.input_file
    sampleid_col = 1
    
    print("CREATE temporary directory")
    temp_dir = tempfile.mkdtemp(prefix="temp_", dir=os.path.abspath(output_maindir))
    logger.info("Temporary files will be written to:{}".format(str(temp_dir)))

    # read mapping file(""sample_id_file", "sample_id_col")
    sample_info = read_mappingfile(input_file, sampleid_col)

    for step in steps_to_run:
        # Step 1: Count number of reads for each sampleids
        if str(step) == '1':
            sample_info = count_reads(sample_info, input_seqdir, output_maindir, step, temp_dir)
        # Step 2: Quality Check of fastq files using fatsqc and multiqc
        elif str(step) == '2':
            sample_info = quality_check(sample_info, input_seqdir, output_maindir, step, temp_dir)
        # Step 3: Run trimmomatic
        elif str(step) == '3':
            logger.debug("Step 3 - Input_seqdir = {}".format(input_seqdir))
            sample_info = trim_reads(sample_info, input_seqdir, output_maindir, step, temp_dir)
        # Step 4: Filter human reads using deconseq
        elif str(step) == '4':
            input_seqdir = os.path.join(output_maindir, "r03_trimmed_reads")
            sample_info = check_output(sample_info, input_seqdir, 3)
            if not os.path.isdir(input_seqdir):
                logger.exception("CRITICAL ERROR: Directory containing trimmomatic ouptut doesn't exists!")
                sys.exit(1)
            else:
                logger.debug("Step 4 - Input_seqdir = {}".format(input_seqdir))
                sample_info = filter_human_reads(sample_info, input_seqdir, output_maindir, step, temp_dir)
        # Step 5: fix paired-end reads
        elif str(step) == '5':
            input_seqdir = os.path.join(output_maindir, "r04_filtered_reads")
            sample_info = check_output(sample_info, input_seqdir, 4)
            if not os.path.isdir(input_seqdir):
                logger.exception("CRITICAL ERROR: Directory containing deconseq ouptut doesn't exists!")
                sys.exit(1)
            else:
                sample_info = fix_paired_reads(sample_info, input_seqdir, output_maindir, step, temp_dir)
        # Step 6: Quality Check of cleaned reads
        elif str(step) == '6':
            input_seqdir = os.path.join(output_maindir, "r05_clean_reads")
            sample_info = check_output(sample_info, input_seqdir, 5)
            if not os.path.isdir(input_seqdir):
                logger.exception("CRITICAL ERROR: Directory containing clean reads doesn't exists!")
                sys.exit(1)
            else:
                sample_info = quality_check(sample_info, input_seqdir, output_maindir, step, temp_dir)
        # Step 7: Count number of reads for each sampleids
        elif str(step) == '7':
            input_seqdir = os.path.join(output_maindir, "r05_clean_reads")
            sample_info = check_output(sample_info, input_seqdir, 5)
            if not os.path.isdir(input_seqdir):
                logger.exception("CRITICAL ERROR: Directory containing clean reads doesn't exists!")
                sys.exit(1)
            else:
                sample_info = count_reads(sample_info, input_seqdir, output_maindir, step, temp_dir)
        else:
            logger.error("Invalid step number! Exiting...", exc_info=True)
            sys.exit(1)

    #Finally: Summarize the results
    generate_summary(sample_info, input_seqdir, output_maindir)

    logger.info("Preprocessing pipeline finished!")


def parse_arguments(argv):

        parser = argparse.ArgumentParser(
            prog = 'preprocess_reads.py',
            formatter_class=argparse.RawTextHelpFormatter,
            description = 'Pipeline for pre-processing of PE sequencing files')

        parser.add_argument(
            '--version',
            action = 'version',
            version = "%(prog)s v" + VERSION)
        parser.add_argument(
            '-i', '--inputfile',
            dest = "input_file",
            required=True,
            help = "Enter input file name that contains list of sample ids\n\n")
        parser.add_argument(
            '-d', '--seqdir',
            dest = "seq_dir",
            required = True,
            help = "Enter path to sequence directory\n\n")
        parser.add_argument(
            '-o', '--outdir',
            dest = "output_dir",
            required = True,
            help = "Enter name of the output dir\n\n")
        parser.add_argument(
            '--step',
            dest = "step_num",
            required = True,
            default = 0,
            type=str,
            help = "Enter the step number separated by ':' to run multiple steps of the pipeline(e.g. 1:6)." \
            " To run all steps, enter 0\n" +
            "[DEFAULT: 0, run all steps in order:]\n\n" +
            "Step 1: Counts the number of reads for each samples\n" +
                    "Step 2: Run fastqc on each sample\n" +
                    "Step 3: Run trimmomatic using used defined parameters. Modify config" +
                        " file 'preprocess_config' to change the parameter settings\n" +
                    "Step 4: Run deconseq to remove human reads.\n" +
                    "Step 5: Fix paired end reads\n" +
                    "Step 6: Run fastqc again on cleaned sequencing reads\n" +
                    "Step 7: Count reads in each sample\n\n")
        parser.add_argument(
            '-c', '--config',
            dest = "config_file",
            required = True,
            help = "Provide the path to config file\n\n")
        parser.add_argument(
            '--log-level',
            default = config.log_level,
            choices=config.log_level_choices,
            help="level of messages to display in log\n" +
            "[DEFAULT: " + config.log_level + "]\n\n")

        return parser.parse_args()

"""
# Check requirements:
    1. Whether required modules are installed and loaded
    2. file formats and permissions
"""
def check_requirements(args):

    # Check that fastqc is installed
    fastqc_path = util.find_exe_in_path("fastqc")
    print "fastqc path = {}".format(str(fastqc_path))
    if not fastqc_path:
                sys.exit("CRITICAL ERROR: The fastqc can not be found. "
                    "Please check whether fastqc module is installed and loaded.")

    # Check that multiqc is installed
    multiqc_path = util.find_exe_in_path("multiqc")
    print "multiqc path = {}".format(str(multiqc_path))
    if not multiqc_path:
                sys.exit("CRITICAL ERROR: The multiqc tool can not be found. "
                    "Please check whether multiqc module is installed and loaded.")

    # Check that trimmomatic is installed
    trim_home = str(os.getenv('TRIMMOMATIC_HOME',""))
    trim_path = util.find_exe_in_path("trimmomatic")
    print("trim_path = {}".format(str(trim_path)))
    if not trim_path:
        sys.exit("CRITICAL ERROR: The trimmomatic can not be found. "
                            "Please check whether trimmomatic module is installed and loaded.")

    #if trim_home:
    #    trim_jar = str(trim_home) + "/trimmomatic-" + str(os.path.basename(trim_home)) + ".jar"
    #    print("trim jar location = {}".format(trim_jar))
    #    proc = subprocess.Popen(['java', '-jar', trim_jar, '-version'], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    #    output = proc.communicate()[0].strip()
    #    print "trimmomatic version = {}".format(str(output))
    #    if not output:
    #        sys.exit("CRITICAL ERROR: The trimmomatic can not be found. "
    #                        "Please check whether trimmomatic module is installed and loaded.")
    #else:
    #    sys.exit("CRITICAL ERROR: The trimmomatic can not be found. "
    #                            "Please check whether trimmomatic module is installed and loaded.")

    # Check that deconseq is installed
    deconseq_path = util.find_exe_in_path("deconseq.pl")
    print "deconseq path = {}".format(str(deconseq_path))
    if not deconseq_path:
        sys.exit("CRITICAL ERROR: The deconseq.pl executable can not be found. "
                    "Please check the install.")


"""
# read config file and write configuration settings based on the arguments
"""
def read_configuration(args):

    args.config_file = os.path.abspath(args.config_file)
    args.input_file = os.path.abspath(args.input_file)

    # Check if config file exists and is readable
    if not os.path.isfile(args.config_file):
        sys.exit("CRITICAL ERROR: Can not find the input file:%s" %(args.config_file))
    if not os.access(args.config_file, os.R_OK):
        sys.exit("CRITICAL ERROR: Not able to read input file selected: " + args.config_file)

    # Check if input file exists and is readable
    if not os.path.isfile(args.input_file):
        sys.exit("CRITICAL ERROR: Can not find the input file:%s" %(args.input_file))
    if not os.access(args.input_file, os.R_OK):
        sys.exit("CRITICAL ERROR: Not able to read input file selected: " + args.input_file)

    # Check the output directory is writable
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create the directory")

    if not os.access(output_dir, os.W_OK):
        sys.exit("CRITICAL ERROR: The output directory is not " +
                    "writeable. This software needs to write files to this directory.\n" +
                    "Please select another directory.")

    print("Output files will be written to: " + output_dir)

    # Configure the logger
    log_file = os.path.join(output_dir , "zz_preprocess.log")
    logging.basicConfig(filename=log_file,format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', level=getattr(logging,args.log_level), filemode='w', datefmt='%m/%d/%Y %I:%M:%S %p')

    #handler = logging.handlers.RotatingFileHandler('zz_analysis.log', maxBytes=2000000, backupCount=20)
    #logger.addHandler(handler)

    # write the version of the software to the log
    logger.info("Running preprocessing pipeline v" + VERSION)
    logger.info("Command: {}".format(" ".join(sys.argv)))

    # Read and write configuration settings
    config_opts = config.read_configfile(args.config_file)
    #config.log_settings()

    return config_opts

"""
# Read mapping file and returns a list.
#
"""
def read_mappingfile(filename, col):

    logger.info("Reading Mapping file: " + filename)
    sample_info = defaultdict(dict)

    try:
        s = os.stat(filename)
        if s.st_size == 0:
            logger.exception("Error Occurred! The file {} is empty".format(filename))
            sys.exit(1)
    except OSError as e:
        logger.exception("Error Occurred: " + str(e))
        sys.exit(2)

    col_index = int(col) - 1

    if col_index < 0:
        logger.error('Invalid column number is specified for read the input file')
        sys.exit(1)

    with open(filename) as fh:
        for line in fh:
            columns = line.rstrip('\n').split('\t')
            if len(columns)> int(col_index):
                sample_info[columns[col_index]] = {1:'N', 2:'N', 3:'N', 4:'N', 5:'N', 6:'N', 7:'N'}
            else:
                logger.error('Invalid column number is specified for read the input file')
                sys.exit(1)

    logger.debug("IDs: {}".format(str(sample_info)))
    #fh.close()

    return sample_info


"""
# This function will count the number of reads in paired end sequencing file
# given the list of sample ids and directory containing sequencing reads
"""
def count_reads(sample_info, input_seqdir, output_maindir, step_num, temp_dir):

    logger.info("Step {}: Counting reads for each sample".format(str(step_num)))

    #missing_ids = []
    f_suffix = config_opts.general['file_suffix_fw']
    r_suffix = config_opts.general['file_suffix_rv']

    if str(step_num) == '1':
    	script_file = os.path.join(output_maindir, "s01_countreads.sh")
	output_file = os.path.join(output_maindir, "r01_readcount.txt")
    elif str(step_num) == '7':
	script_file = os.path.join(output_maindir, "s07_countreads.sh")
	output_file = os.path.join(output_maindir, "r07_readcount.txt")
    else:
	logger.exception("count_reads - Invalid step number!")
	sys.exit(1)	
 
    ## check if a file exists##
    ## if exists, delete it else show message on screen ##
    if os.path.exists(output_file):
        try:
            os.remove(output_file)
        except OSError, e:
            logger.exception("Error: %s - %s." % (e.output_file,e.strerror))
   
    fh1 = open(script_file, 'w')

    for id, data in sample_info.iteritems():

	logger.debug("data = {}".format(str(data)))
    	if str(step_num) == '1':
            fw_readfile = os.path.join(input_seqdir, "{}_{}.fastq".format(str(id), str(f_suffix)))
            rv_readfile = os.path.join(input_seqdir, "{}_{}.fastq".format(str(id), str(r_suffix)))

	elif str(step_num) == '7':
            fw_readfile = os.path.join(input_seqdir, "{}_{}_CLEAN.fastq".format(str(id), str(f_suffix)))
            rv_readfile = os.path.join(input_seqdir, "{}_{}_CLEAN.fastq".format(str(id), str(r_suffix)))

	else:
            logger.exception("Invalid step number")
            sys.exit(1)

    	if os.access(fw_readfile, os.R_OK) and os.access(rv_readfile, os.R_OK):
	    logger.debug("Counting number of reads in {}".format(fw_readfile))
            #cmd = 'fgrep -h -c "`head -n 1 ' + fw_readfile + ' | awk -F\'[:]\' \'{print $1}\'`" ' \
            #    + str(os.path.join(input_seqdir, str(id) + "_*.fastq" )) \
            #    + ' | xargs echo ' + id + ' | sed \'s/ /\\t/g\' >> ' + output_file + '; echo exit_code:$?'
            cmd = 'cnts=`echo $(cat {} | wc -l)/4|bc`; echo -e "{}\t$cnts" >> {}; echo exit_code:$?'.format(str(fw_readfile), str(id), str(output_file))
            logger.debug("cmd = {}".format(str(cmd)))
            #cnts=$(echo $(cat $i | wc -l)/4 | bc); echo -e "${i}\t${cnts}"
            sample_info[id][1] = 'Y'
            fh1.write(cmd + "\n")

    	else:
            logger.exception("Either file {} is missing or is not readable".format(rv_readfile))
            sample_info[id][1]="N"
            continue


    fh1.close()

    if os.stat(script_file).st_size > 0:
        sbatch_params = '--mem=16000M'
        job_status = util.run_sbatch_script(script_file, 10, output_maindir, sbatch_params, temp_dir)
        logger.info("job status = {}".format(str(job_status)))
        #sample_info = check_output(sample_info, output_subdir, step_num)
        #pass

    else:
        logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
        sys.exit(1)

    #if len(missing_ids) > 0:
    #   logger.info("missing_ids = " + str(missing_ids))

    logger.debug("sample_info = {}".format(str(sample_info)))
    return sample_info


def quality_check(sample_info, input_seqdir, output_maindir, step_num, temp_dir):

    logger.info("Step: " + str(step_num) + " Quality check using fastqc")

    output_subdir = os.path.join(output_maindir, "r0" + str(step_num) + "_fastqc_output")
    script_file = os.path.join(output_maindir, "s0" + str(step_num) + "-1_run_fastqc.sh")

    if not os.path.isdir(output_subdir):
                try:
                        os.mkdir(output_subdir)
                except EnvironmentError:
                        sys.exit("CRITICAL ERROR: Unable to create the directory")

    fh2 = open(script_file,'w')

    index = int(step_num) - 1

    for id, data in sample_info.iteritems():

        #if data[index] == 'Y':
        cmd = "fastqc -o %s %s" %(str(output_subdir), str(os.path.join(input_seqdir, id + "_*.fastq")))
        logger.debug("cmd = {}".format(str(cmd)))
        fh2.write(cmd + '\n')
        #else:
        #   sample_info[id][step_num] = 'N'

    fh2.close()

    if os.stat(script_file).st_size > 0:
        sbatch_params = '--mem=16000M '
        job_status = util.run_sbatch_script(script_file, 1, output_maindir, sbatch_params, temp_dir)
        logger.info("job status = {}".format(str(job_status)))

        if job_status == 'SUCCESS':
            cmd2 = "multiqc --interactive -m fastqc -f -o {} {}".format(str(os.path.join(output_subdir, "zz_multiqc")), str(output_subdir))
            logger.debug("cmd2 = {}".format(str(cmd2)))	
    
            script_file2 = os.path.join(output_maindir, "s0" + str(step_num) + "-2_run_multiqc.sh")
            fh3 = open(script_file2, 'w')
            fh3.write(cmd2 + '\n')
            fh3.close()
            sbatch_params = '--mem=16000M '
            cmd_status = util.run_sbatch_script(script_file2, 1, output_maindir, sbatch_params, temp_dir)
    
    		# cmd_status = subprocess.check_call(cmd)
            logger.debug("cmd_status = {}".format(str(cmd_status)))
            sample_info = check_output(sample_info, output_subdir, step_num)
        else:
            logger.exception('Error occurred at quality check step')
            sys.exit(1)

    else:
        logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
        sys.exit(1)

    logger.debug("sample info = {}".format(str(sample_info)))

    return sample_info



def trim_reads(sample_info, input_seqdir, output_maindir, step_num, temp_dir):

    logger.info("Step 3: Run Trimmomatic")

    f_suffix = config_opts.general['file_suffix_fw']
    r_suffix = config_opts.general['file_suffix_rv']
    output_subdir = os.path.join(output_maindir,"r03_trimmed_reads")
    script_file = os.path.join(output_maindir,"s03_run_trimmomatic.sh")

    #trim_home = str(os.getenv('TRIMMOMATIC_HOME',""))
    trim_exe = str(util.find_exe_in_path("trimmomatic"))
    if trim_exe:
        trim_path =  os.path.dirname(os.path.realpath(trim_exe))
        trim_home = os.path.dirname(os.path.realpath(trim_path))
        #jardir = os.path.join(trim_home, 'trimmomatic-' + str(os.path.basename(trim_home)) + ".jar")
        logger.debug("trimmomatic home dir = {}".format(trim_home))
    else:
        logger.error("Trimmomatic not installated or in the PATH environemnt")
        sys.exit(1)

    if config_opts.trimmomatic['adapt']:
        logger.debug('trimmomatic adapter information = {}'.format(config_opts.trimmomatic['adapt']))
        if config_opts.trimmomatic['adapt'] == "nextera":
            adapt_path = str(os.path.join(trim_home, "share", "adapters", 'NexteraPE-PE.fa'))
        elif config_opts.trimmomatic['adapt'] == "truseq":
            adapt_path = str(os.path.join(trim_home, "share", "adapters", 'TruSeq3-PE-2.fa'))
        else:
            logger.error("Invalid trimmommatic adapter information provided! Exiting..")
            sys.exit(1)
    else:
        logger.error("Invalid trimmommatic adapter information provided! \
            Modify the adapter information under trimmomatic section in td_config.cfg file. \
            Available options = nextera|truseq")
        sys.exit(1)

    logger.debug("trimmomatic adapt = {}".format(str(adapt_path)))

    leading = config_opts.trimmomatic['leading']
    trailing = config_opts.trimmomatic['trailing']
    illuminaclip = config_opts.trimmomatic['illuminaclip']
    slidingwindow = config_opts.trimmomatic['sliding_window']
    minlen = config_opts.trimmomatic['minlen']

    if not os.path.isdir(output_subdir):
        try:
            os.mkdir(output_subdir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create the directory")

    fh3 = open(script_file,'w')

    logger.debug("sample_info = ".format(str(sample_info)))

    for id, data in sample_info.iteritems():

    	logger.debug("data = {}".format(str(data)))
        #if data[2] == 'Y':
    	fwfile = os.path.join(input_seqdir, "{}_{}.fastq".format(str(id), str(f_suffix)))
    	rvfile = os.path.join(input_seqdir, "{}_{}.fastq".format(str(id), str(r_suffix)))
    	paired_fw = os.path.join(output_subdir, "{}_{}_TRIMMED_PAIRED.fastq".format(str(id), str(f_suffix)))
    	paired_rv = os.path.join(output_subdir, "{}_{}_TRIMMED_PAIRED.fastq".format(str(id), str(r_suffix)))
    	unpaired_fw = os.path.join(output_subdir, "{}_{}_TRIMMED_UNPAIRED.fastq".format(str(id), str(f_suffix)))
    	unpaired_rv = os.path.join(output_subdir, "{}_{}_TRIMMED_UNPAIRED.fastq".format(str(id), str(r_suffix)))


        #this command has been changed to include a call to the amount of memory required. (Set to 4G)
        cmd = "JAVA_ARGS=\"-Xmx4096m\" trimmomatic" + " PE -phred33 " + fwfile + " " \
            + rvfile + " " + paired_fw + " " + unpaired_fw + " " + paired_rv + " " + unpaired_rv \
            + " ILLUMINACLIP:"+ str(adapt_path) + ":" + str(illuminaclip) + " LEADING:" + str(leading) \
            + " TRAILING:" + str(trailing) + " SLIDINGWINDOW:" + str(slidingwindow) + " MINLEN:" + str(minlen)

        #cmd = "trimmomatic" + " PE -phred33 " + fwfile + " " \
        #    + rvfile + " " + paired_fw + " " + unpaired_fw + " " + paired_rv + " " + unpaired_rv \
        #    + " ILLUMINACLIP:"+ str(adapt_path) + ":" + str(illuminaclip) + " LEADING:" + str(leading) \
        #    + " TRAILING:" + str(trailing) + " SLIDINGWINDOW:" + str(slidingwindow) + " MINLEN:" + str(minlen)

    	#cmd = "java -Xms1024m -Xmx1024m -jar " + jardir + " PE -phred33 " + fwfile + " " \
        #    + rvfile + " " + paired_fw + " " + unpaired_fw + " " + paired_rv + " " + unpaired_rv \
        #    + " ILLUMINACLIP:"+ str(adapt_path) + ":" + str(illuminaclip) + " LEADING:" + str(leading) \
        #    + " TRAILING:" + str(trailing) + " SLIDINGWINDOW:" + str(slidingwindow) + " MINLEN:" + str(minlen)

    	logger.debug("cmd = {}".format(str(cmd)))

    	fh3.write(cmd + '\n')

    fh3.close()

    if os.stat(script_file).st_size > 0:
        sbatch_params = '--mem=16000M'
        util.run_sbatch_script(script_file, 1, output_maindir, sbatch_params, temp_dir)
        sample_info = check_output(sample_info, output_subdir, step_num)
        #pass

    else:
        logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
        sys.exit(1)

    logger.debug("sample info = {}".format(str(sample_info)))

    return sample_info


"""
# Filter out human contamination using deconseq
#
"""
def filter_human_reads(sample_info, input_seqdir, output_maindir, step_num, temp_dir):

    logger.info("Step 4: Filter human reads using deconseq")
    f_suffix = config_opts.general['file_suffix_fw']
    r_suffix = config_opts.general['file_suffix_rv']
    output_subdir = os.path.join(output_maindir,"r04_filtered_reads")
    script_file = os.path.join(output_maindir,"s04_run_deconseq.sh")
    #deconseq_path = util.find_exe_in_path("deconseq.pl")
    deconseq_path = config_opts.deconseq['deconseq_path']
    if not deconseq_path:
        logger.warning("Deconseq is not loaded in the module. \
                We will use the default option from config file")
        deconseq_path = config_opts.deconseq['deconseq_path']

    deconseq_db = config_opts.deconseq['deconseq_db']
    #suffix = deconseq_out_suffix

    if not os.path.isdir(output_subdir):
        try:
            os.mkdir(output_subdir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create the directory")

    fh4 = open(script_file,'w')

    for id, data in sample_info.iteritems():
    	if data[3] == 'Y':

            fwin = os.path.join(input_seqdir, "{}_{}_TRIMMED_PAIRED.fastq".format(str(id), str(f_suffix)))
            rvin = os.path.join(input_seqdir, "{}_{}_TRIMMED_PAIRED.fastq".format(str(id), str(r_suffix)))
            fwout= os.path.join(output_subdir, "{}_{}_TRIMMED_PAIRED_FILTERED".format(str(id), str(f_suffix)))
            rvout= os.path.join(output_subdir, "{}_{}_TRIMMED_PAIRED_FILTERED".format(str(id), str(r_suffix)))
            cmd = "{} -f {} -o {} -dbs {} -id {} ; {} -f {} -o {} -dbs {} -id {}" \
                    .format(deconseq_path, fwin, fwout, deconseq_db, str(id), deconseq_path, rvin, rvout, deconseq_db, str(id))

            logger.debug("cmd = {}".format(str(cmd)))
                        #sample_info[id][4] = 'Y'

            fh4.write(cmd + '\n')
        else:
            sample_info[id][4] = 'N'

    fh4.close()

    if os.stat(script_file).st_size > 0:
        sbatch_params = '--mem=60000M'
        util.run_sbatch_script(script_file, 1, output_maindir, sbatch_params, temp_dir)
        sample_info = check_output(sample_info, output_subdir, step_num)
        #pass

    else:
        logger.error("The script file {} is empty. Exiting....".format(script_file),exc_info=True)
        sys.exit(1)

    logger.debug("sample info = {}".format(str(sample_info)))

    return sample_info


"""
# Fix paired end reads using Molly's script
#
"""
def fix_paired_reads(sample_info, input_seqdir, output_maindir, step_num, temp_dir):

    logger.info("Step 5: Fix paired end reads")
    f_suffix = config_opts.general['file_suffix_fw']
    r_suffix = config_opts.general['file_suffix_rv']
    output_subdir = os.path.join(output_maindir,"r05_clean_reads")
    script_file = os.path.join(output_maindir,"s05_fix_pairedend.sh")

    if not os.path.isdir(output_subdir):
        try:
            os.mkdir(output_subdir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create the directory")

    fh5 = open(script_file,'w')

    logger.debug("Step 5 sample info =>{}".format(str(sample_info)))

    for id, data in sample_info.iteritems():
    	if data[4] == 'Y':
            fwin = os.path.join(input_seqdir, "{}_{}_TRIMMED_PAIRED_FILTERED".format(str(id), str(f_suffix)), str(id) + "_clean.fq")
            rvin = os.path.join(input_seqdir, "{}_{}_TRIMMED_PAIRED_FILTERED".format(str(id), str(r_suffix)), str(id) + "_clean.fq")
            fwout = os.path.join(output_subdir,"{}_{}_CLEAN.fastq".format(str(id), str(f_suffix)))
            rvout = os.path.join(output_subdir, "{}_{}_CLEAN.fastq".format(str(id), str(r_suffix)))
	    singletons = os.path.join(output_subdir, "{}_SINGLETON_CLEAN.fastq".format(str(id)))

            cmd = "repair.sh ow=t in={in_fw} in2={in_rv} out={out_fw} out2={out_rv} outs={out_sl} repair=t"\
            .format(
                in_fw = fwin,
                in_rv = rvin,
                out_fw = fwout,
                out_rv = rvout,
                out_sl = singletons
            )
            #cmd = "fix_paired_end_MB.py " + fwin + " " + rvin + " " + fwout + " " + rvout
            logger.debug("cmd = {}".format(str(cmd)))

            sample_info[id][5] = 'Y'
            fh5.write(cmd + '\n')
        else:
            sample_info[id][5] = 'N'

    fh5.close()

    if os.stat(script_file).st_size > 0:
        sbatch_params = '--mem=30000M'
        util.run_sbatch_script(script_file, 1, output_maindir,sbatch_params, temp_dir)
        sample_info = check_output(sample_info, output_subdir, step_num)
        #pass

    else:
        logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
        sys.exit(1)

    logger.debug("sample info = {}".format(str(sample_info)))

    return sample_info


"""
# Check for output files after each step
#
"""
def check_output(sample_info, output_dir, step_num):

    logger.debug("Checking output - step {}".format(str(step_num)))
    logger.debug("sample info ->{}".format(str(sample_info)))
    f_suffix = config_opts.general['file_suffix_fw']
    r_suffix = config_opts.general['file_suffix_rv']
    suffix_ext = {1:'.fastq', 2:'fastqc.html', 3:'TRIMMED_PAIRED.fastq' , 4:'TRIMMED_PAIRED_FILTERED' , 5:'CLEAN.fastq' , 6:'CLEAN_fastqc.html'}
    step_num = int(step_num)

    for id in sample_info.keys():
        fw_file = os.path.join(output_dir, '{}_{}_{}'.format(str(id), str(f_suffix), str(suffix_ext[int(step_num)])))
        rv_file = os.path.join(output_dir, '{}_{}_{}'.format(str(id), str(r_suffix), str(suffix_ext[int(step_num)])))

        logger.debug("Checking output - fw_file = {} and rv_file = {}".format(fw_file, rv_file))

        if os.access(fw_file, os.R_OK) and os.access(rv_file, os.R_OK):
            sample_info[id][step_num] = 'Y'

        else:
            logger.exception("Either forward or reverse file for {} is missing or is not readable".format(id))
            sample_info[id][step_num] = 'N'


    return sample_info


"""
# Summarize results
"""
def generate_summary(sample_info, input_seqdir, output_maindir):

    logger.info("Summarizing results")
    output_file = os.path.join(output_maindir, 'r08_summary.txt')

    fh7 = open(output_file, 'w')
    fh7.write("#SampleID\tStep_1\tStep_2\tStep_3\tStep_4\tStep_5\tStep_6\tStep_7\n")

    for id, data in sample_info.iteritems():
        fh7.write(id + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] \
                + "\t" + data[5] + "\t" + data[6] + "\t" + data[7] + "\n")


if __name__ == "__main__":
        main(sys.argv[1:])
