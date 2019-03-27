#!/usr/bin/env python

"""
File Name    : sequence_assembly_MB.py
Author       : Manish Boolchandani, manish@wustl.edu
Created On   : 2017-04-24
Last Modified: 2017-04-25

Description  :	Performs assembly using Spades, IDBA-UD and Velvet 
		and also runs the quast to evaluate genome assembly
		Final step : Generate summary for each samples

Dependencies: 	Spades, Velvet, IDBA_UD, QUAST

Usage: sequence_assembly.py [-h] [--version] -i INPUT_FILE -d SEQ_DIR -o OUTPUT_DIR --step STEP_NUM
                           [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
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
from collections import defaultdict

# Other imports
import main_config_MB as config
import utilities_MB as util

#Version
VERSION = "0.1.0"

#Logger configuration and setup
logger = logging.getLogger(__name__)

# Maximum steps in the pipeline

assembly_list = []

def main(argv):

	logger.info("Started Preprocessing Pipeline v" + VERSION)
	
	# Parse arguments from command line
	args = parse_arguments(argv)
	#print(args)

	if not 'all' in args.assembly_list:
		assembly_list = args.assembly_list
	else:
		assembly_list = ['velvet','spades','idba']

	# Check the required files, software and dependencies  
	#check_requirements(args)	
	
	# Update and write the config settings to the log file
	update_configuration(args)
	
	# Write configuration settings
	config.log_settings()
		
	
	logger.info("assembler list = {}".format(str(assembly_list)))
 
	output_maindir = args.output_dir
	input_seqdir = args.seq_dir
	input_file = args.input_file
	sample_info = read_mappingfile(input_file,3)
	
	for assembler in assembly_list:
		
		if assembler == 'spades':		#Run SPAdes assember
			logger.info("Will now perform assembly using Spades")
			sample_info = run_spades(sample_info, input_seqdir, output_maindir)
		
		elif assembler == 'velvet':		# Perform Velvet assembly using VelvetOptimiser
			logger.info("Will now perform assembly using Velvet")
			sample_info = run_velvet(sample_info, input_seqdir, output_maindir)


		elif assembler == 'idba':		# Perform assembly using IDBA-UD
			logger.info("Will now perform assmebly using idba")
			#logger.debug("Step 3 - Input_seqdir = {}".format(input_seqdir))
			#sample_info = run_idba(sample_info, input_seqdir, output_maindir, step)
				
		else:
	
			logger.error("Invalid assembler name! Exiting...", exc_info=True)
			sys.exit(1) 
	
	##Finally: Summarize the results
	#generate_summary(sample_info, input_seqdir, output_maindir)
	
	logger.info("Assembly pipeline finished!")			


def parse_arguments(argv):
	
	parser = argparse.ArgumentParser(
		prog = 'preprocess_reads.py', 
		description = 'Pipeline for preprocessing of  PE fastq files')
        parser.add_argument(
		'--version', 
		action = 'version', 
		version = "%(prog)s v" + VERSION)
        parser.add_argument(
		'-i', '--inpufile', 
		dest = "input_file", 
		required=True,
		help = "Enter input file name that contains list of sample ids")
        parser.add_argument(
		'-d', '--seqdir', 
		dest = "seq_dir", 
		required = True,
		help = "path to sequence directory")
        parser.add_argument(
		'-o', '--outdir', 
		dest = "output_dir", 
		required = True,
		help = "name of the output dir")
	parser.add_argument(
		'-a', '--assembly',
		dest = "assembly_list",
		required = True,
		choices = ['velvet','spades','idba','all'],
		action='append',
		#nargs = '+',
		help = "Enter the step number to run in the pipeline")
	parser.add_argument(
		'--log-level',
		default = config.log_level,
		choices=config.log_level_choices,
		help="level of messages to display in log\n" +
		"[DEFAULT: " + config.log_level + "]") 

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
	
	# Check that trimmomatic is installed
	trim_home = str(os.getenv('TRIMMOMATIC_HOME',""))
	if trim_home:
		trim_jar = str(trim_home) + "/trimmomatic-" + str(os.path.basename(trim_home)) + ".jar"
        	proc = subprocess.Popen(['java', '-jar', trim_jar, '-version'], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
        	output = proc.communicate()[0].strip()	
		print "trimmomatic version = {}".format(str(output))
		if not output:
                	sys.exit("CRITICAL ERROR: The trimmomatic can not be found. "
                    		"Please check whether trimmomatic module is installed and loaded.")
	else:
		sys.exit("CRITICAL ERROR: The trimmomatic can not be found. "
                                "Please check whether trimmomatic module is installed and loaded.")

	# Check that deconseq is installed
	deconseq_path = util.find_exe_in_path("deconseq.pl")
	print "deconseq path = {}".format(str(deconseq_path)) 		
	if not deconseq_path:
                sys.exit("CRITICAL ERROR: The deconseq.pl executable can not be found. "
                    "Please check the install.")

	# Check that the metaphlan2 executable can be found
        if not util.find_exe_in_path("metaphlan2.py"):
                sys.exit("CRITICAL ERROR: The metaphlan2.py executable can not be found. "
                    "Please check the install.")



"""	
# update the configuration settings based on the arguments
"""
def update_configuration(args):

	args.input_file = os.path.abspath(args.input_file)
	
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
	log_file = os.path.join(output_dir , "zz_assembly.log")

	logging.basicConfig(filename=log_file,format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', level=getattr(logging,args.log_level), filemode='w', datefmt='%m/%d/%Y %I:%M:%S %p')

	#handler = logging.handlers.RotatingFileHandler('zz_analysis.log', maxBytes=2000000, backupCount=20)
	#logger.addHandler(handler)
		
	# write the version of the software to the log
	logger.info("Running sequence assenbly pipeline v" + VERSION)
	logger.info("Command: {}".format(" ".join(sys.argv)))


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
			if not line.startswith('#'):
				columns = line.split('\t')
				if len(columns)> int(col_index):
					sample_info[columns[col_index]] = {1:'N', 2:'N', 3:'N', 4:'N'}
				else:
					logger.error('Invalid column number is specified for read the input file')
					sys.exit(1)
	
			
	logger.info("IDs: {}".format(str(sample_info)))

        return sample_info


"""
# This function will perform velvet assembly on paired end sequencing file
# The funtion takes the list of sample ids and directory containing sequencing reads
# as an input
"""
def run_velvet(sample_info, input_seqdir, output_maindir):
	
	logger.info("Step: Running velvet assembly")
	
	output_subdir = os.path.join(output_maindir,"r02_velvetAssembly")
	if not os.path.isdir(output_subdir):
                try:
                        os.mkdir(output_subdir)
                except EnvironmentError:
                        sys.exit("CRITICAL ERROR: Unable to create the directory")

	script_file = os.path.join(output_maindir, "s02_velvet_script.sh")
	fh1 = open(script_file, 'w')

	velvet_kmer_start = config.velvet_kmer_start
        velvet_kmer_end = config.velvet_kmer_end
        velvet_kmer_step = config.velvet_kmer_step

	velvet_opts = '-t 8 --verbose'

	for id, data in sample_info.iteritems():

		logger.debug("data = {}".format(str(data)))
		fw_readfile = os.path.join(input_seqdir, str(id) + "_FW_CLEAN.fastq")
		rv_readfile = os.path.join(input_seqdir, str(id) + "_RV_CLEAN.fastq")
		
		if os.access(fw_readfile, os.R_OK) and os.access(rv_readfile, os.R_OK):
				
        		velvet_dir_final = os.path.join(output_subdir,str(id))
			velvet_prefix = str(config.velvet_prefix + "_" + str(id))
			
			cmd = "VelvetOptimiser.pl -s {} -e {} -f \'-fastq -shortPaired -separate {} {}\' -t 8 --verbose -d {} -p {}" \
                                .format(velvet_kmer_start, velvet_kmer_end, fw_readfile, rv_readfile, velvet_dir_final, velvet_prefix)
			logger.debug("cmd = {}".format(str(cmd)))
			
			sample_info[id][1] = 'Y'
			fh1.write(cmd + "\n")
 
		else:
	        	logger.exception("Either file {} is missing or is not readable".format(rv_readfile))
			sample_info[id][1]="N"
	                continue

			
	fh1.close() 
		
	if os.stat(script_file).st_size > 0:
		sbatch_params = '--cpus-per-task=8 --mem=30000M'
		util.run_sbatch_script(script_file, 1, output_maindir, sbatch_params)
		sample_info = check_output(sample_info, output_subdir, 'velvet')
		#pass
		
	else:
		logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
		sys.exit(1)
	
	logger.info("sample_info = {}".format(str(sample_info)))	
	return sample_info


"""
# This function will perform assembly using SPAdes
# The funtion takes the list of sample ids and directory containing sequencing reads
# as an input
"""
def run_spades(sample_info, input_seqdir, output_maindir):
	
	logger.info("Step: Running Spades assembly")
	
	output_subdir = os.path.join(output_maindir,"r01_spadesAssembly")
	if not os.path.isdir(output_subdir):
                try:
                        os.mkdir(output_subdir)
                except EnvironmentError:
                        sys.exit("CRITICAL ERROR: Unable to create the directory")

	script_file = os.path.join(output_maindir, "s01_spades_script.sh")
	fh1 = open(script_file, 'w')

	spades_kmers = config.spades_kmers
        spades_opts = config.spades_opts

	for id, data in sample_info.iteritems():

		logger.debug("data = {}".format(str(data)))
		fw_readfile = os.path.join(input_seqdir, str(id) + "_FW_CLEAN.fastq")
		rv_readfile = os.path.join(input_seqdir, str(id) + "_RV_CLEAN.fastq")
		
		if os.access(fw_readfile, os.R_OK) and os.access(rv_readfile, os.R_OK):
							
        		spades_dir_final = os.path.join(output_subdir,str(id))
			cmd = "spades.py -k {} {} --pe2-1 {} --pe2-2 {} -o {}" \
                                .format(spades_kmers, spades_opts, fw_readfile, rv_readfile, spades_dir_final)
			
			logger.debug("cmd = {}".format(str(cmd)))
			
			sample_info[id][1] = 'Y'
			fh1.write(cmd + "\n")
 
		else:
	        	logger.exception("Either file {} is missing or is not readable".format(rv_readfile))
			sample_info[id][1]="N"
	                continue

			
	fh1.close() 
				
	if os.stat(script_file).st_size > 0:
		sbatch_params = '--mem=30000M --cpus-per-task=8'
		util.run_sbatch_script(script_file, 1, output_maindir, sbatch_params)
		sample_info = check_output(sample_info, output_subdir, 'spades')
		#pass
		
	else:
		logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
		sys.exit(1)
	
	logger.info("sample_info = {}".format(str(sample_info)))	
	return sample_info

"""
# Check for output files after each step
#
"""
def check_output(sample_info, output_dir, assembler):

	logger.debug("Checking output - {}".format(str(assembler)))

	suffix = {'spades':'contigs.fasta', 'velvet':'contigs.fa', 'idba':'contigs.fasta' , 'abyss':'contigs.fasta'}  
	
	for id in sample_info.keys():
		assembly_file = os.path.join(output_dir, id, str(suffix[assembler]))
		
		logger.info("Checking output - {}".format(assembly_file))
		
		if os.access(assembly_file, os.R_OK):
                	sample_info[id][assembler] = 'Y'

      		else:
                	logger.exception("Either forward or reverse file for {} is missing or is not readable".format(id))
                	sample_info[id][assembler] = 'N'


	return sample_info


"""
# Summarize results
"""
def generate_summary(sample_info, input_seqdir, output_maindir):

	logger.info("Summarizing results")
	output_file = os.path.join(output_maindir, 'r07_summary.txt')
	
	fh7 = open(output_file, 'w')
	fh7.write("#SampleID\tSPAdes\tVelvet\tIDBA\tAbyss\n")
	
	for id, data in sample_info.iteritems():
		fh7.write(id + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] + "\n") 


if __name__ == "__main__":
        main(sys.argv[1:])
