#!/usr/bin/env python

"""
File Name    : subsample_PEreads_MB.py
Author       : Manish Boolchandani, manish@wustl.edu
Created On   : 2016-07-05
Last Modified: 2016-07-05

Description  :  This program will take the list of sample ids, path to 
                sequence file (raw reads) and will perform subsampling 
                based on number of reads provided by the user.

Dependencies:   

Usage:  subsample_PEreads.py [-h] [-v] [-i INPUTFILE] [-d SEQDIR]
                            [-m MIN_READS] [-x MAX_READS] [-s STEP_NUM]
                            [-o OUTPUTDIR]
CHANGE LOG:
TODO:

"""

# Python imports
import os, operator, sys, time, re, logging, argparse, subprocess, glob, random
import logging.handlers
from collections import defaultdict

#Logger configuration and setup
LEVEL = logging.DEBUG
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger(__file__)
logger.setLevel(LEVEL)
handler = logging.handlers.RotatingFileHandler('y2_subsample.log', maxBytes=2000000, backupCount=20)
handler.setLevel(LEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)

#Version
VERSION = "1.0.0"
##########################################################################################
# Main Program
#
##########################################################################################

def main(argv):

	logger.info("STARTED PIPELINE!!")
	
	parser = argparse.ArgumentParser(prog = 'subsample_PEreads.py', description = 'A program to subsample PE reads ')
	parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s v' + VERSION)
	parser.add_argument('-i', '--ifile', dest = "inputfile", help = "Enter input file name that contains list of sample ids")
	parser.add_argument('-d', '--seqdir', dest = "seqdir", help = "path to sequence directory")
	parser.add_argument('-m', '--min', dest = "min_reads", help = "The minimum number of reads to sart with")
	parser.add_argument('-x', '--max', dest = "max_reads", help = "The maximum number of reads to samplek")
	parser.add_argument('-s', '--step', dest = "step_num", default = 0, help = "Please enter step size to increment from the lowest to highest read depth")
	parser.add_argument('-o', '--outdir', dest = "outputdir", help = "name of the output dir")

	args = parser.parse_args()

	# check to see that input file is provided
	if not args.inputfile:
		logger.exception("You must provide input file containing sample ids. Check usage with subsample_PEreads.py -h")
		parser.exit(status=1, message="You must provide input file containing sample ids\nCheck usage with 'subsample_PEreads.py -h'.\n\n")
	else:
		inputfile = args.inputfile

	if not args.outputdir:
		outputdirname = os.path.splitext(os.path.basename(inputfile))[0] + "_out"
		outdir = os.path.join(os.getcwd(), outputdirname)
		logger.warning("Output directory is not provided. Writing the result to {} " . format(outputdir))
	else:
		outdir = os.path.realpath(args.outputdir)

	if not args.seqdir:
		logger.exception("You must provide path to directory that contains sequences or metaphlan output")
		parser.exit(status=1, message="You must provide path to sequence directory\n")
	else:
		seqdir = os.path.realpath(args.seqdir) 

	if not args.min_reads:
		min_reads = 100000
	else:
		min_reads = args.min_reads
	
	if not args.max_reads:
		max_reads = 1000000
	else:
		max_reads = args.max_reads
	
	if not args.step_num:
		step_num = 100000
	else:
		step_num = args.step_num
	
	logger.debug("minimum read depth: " + min_reads)
	logger.debug("maximum read depth: " + max_reads)
	logger.debug("step_num: " + step_num)
	
	try:
		
		if inputfile:

			metaphlan_output = ""
			logger.info("Input filename = " + inputfile) 

			if not os.path.isdir(outdir): subprocess.check_call(['mkdir', outdir]);
				
			missing_ids = []
			idList = []
			#logger.info("STEP 1: Run Metaphlan")

			idList = read_mappingfile(inputfile)
			#logger.debug("data in input file {}".format(str(dict)))
                        scriptfile = os.path.join(outdir, "slurm_subsample.sh")
		        fh1 = open(scriptfile, 'w')
			mappingfile = os.path.join(outdir, "e01_mappingfile.txt")
			fh2 = open(mappingfile,'w')
                        
			for id in idList:
				#logger.debug("id:" + str(id))
				fw_readfile = os.path.join(seqdir, str(id) + "_FW_CLEAN.fastq")
				rv_readfile = os.path.join(seqdir, str(id) + "_RV_CLEAN.fastq")
				
				if os.path.isfile(fw_readfile) and os.access(fw_readfile, os.R_OK):
                                        #logger.debug("File exists and is readable")
                                        pass
                                else:
                                        logger.exception("Either file {} is missing or is not readable".format(fw_readfile))
                                        missing_ids.append(os.path.splitext(os.path.basename(fw_readfile))[0])
                                        continue

                                if os.path.isfile(rv_readfile) and os.access(rv_readfile, os.R_OK):
                                        #logger.debug("File exists and is readable")
                                        pass
                                else:
                                        logger.exception("Either file {} is missing or is not readable".format(rv_readfile))
                                        missing_ids.append(os.path.splitext(os.path.basename(rv_readfile))[0])
                                        continue
				
				logger.info("Counting number of reads in {}".format(fw_readfile))
				cmd = 'LC_ALL=C fgrep -c "`head -n 1 ' + fw_readfile + " | awk -F '[:]' '{print $1}'`\" " + fw_readfile
				#logger.debug("cmd = "+ cmd)
				proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
        			fwreadcnt = int(proc.communicate()[0].strip())

				cmd = 'LC_ALL=C fgrep -c "`head -n 1 ' + rv_readfile + " | awk -F '[:]' '{print $1}'`\" " + rv_readfile
				#logger.debug("cmd = "+ cmd)
                                proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
                                rvreadcnt = int(proc.communicate()[0].strip())

        			#output = subprocess.check_output(cmd, shell=True)
        			if fwreadcnt == rvreadcnt:
					total_reads = fwreadcnt
					data = subsample_reads(id, fw_readfile, rv_readfile, total_reads, min_reads, max_reads, step_num, outdir)
				 	for key , value  in data.items():
						fh1.write(value['cmd'] + "\n")
						fh2.write("{}\t{}\t{}\t{}\t{}\n".format( id + "_" + str(key), id, key, value['readcnt'], total_reads))
				else:
					logger.exception("forward read count is not equal to reverse read count")

			fh1.close() 
			fh2.close()

			if os.stat(scriptfile).st_size > 0:
		                srun_script(scriptfile)
        		else:
                		logger.exception("The file {} is empty. Exiting....".format(scriptfile))
                		sys.exit()
	
			if len(missing_ids) > 0:
				logger.info("missing_ids = " + str(missing_ids))

		else:
			logger.exception("Please input mapping file containing sampleIds")
			sys.exit(1)

	except NameError as e:
		print "Error Occurred:", e
		logger.exception("Error Occurred: " + str(e))
		sys.exit(1) 
	except:
		logger.exception("Error Occurred")
		sys.exit(1)

	logger.info("THE END!!\n")
	
	return;


#################################################################################################
# 
# Read mapping file and return a dictionary. 
# Ex. for each row store [id][SeqsampleId]=[ReadCount]
#
#################################################################################################
def read_mappingfile(filename):

        logger.debug("Reading Mapping file: " + filename)
	#logger.debug("Step Num: " + str(step_num))

        try:
                s = os.stat(filename)
                if s.st_size == 0:
                        #print "The file {} is empty".format(filename)
                        logger.exception("Error Occurred! The file {} is empty".format(filename))
                        sys.exit(1)
        except OSError as e:
                #print e
                logger.exception("Error Occurred: " + str(e))
                sys.exit(2)

        fh = open(filename, "r")
        first_line = fh.readline()

        sampleInfo = [line.rstrip('\n').split("\t") for line in fh]
	idList = []
	for id in sampleInfo:
		idList.extend(id)
	
	logger.info("Total IDs: {}".format(len(idList)))
        
	fh.close()

        return idList

#########################################################################################
#
# Subsample Paired end reads
#
#
#########################################################################################

def subsample_reads(id, fw_readfile, rv_readfile, total_reads, min_reads, max_reads, step_num, outdir):
	sampleid = id
	total_reads = int(total_reads)
	min = int(min_reads)
	max = int(max_reads)
	step = int(step_num)
        cmdList = []
	
	logger.debug("Number of reads for " + sampleid + " = " + str(total_reads) )
	
	#scriptfile = os.path.join(outdir, "s01_subsample.sh")
	#fh = open(scriptfile, 'w')
	
	if max > total_reads:
		max = total_reads

	data = {}
	cnt = 1;
	for i in range(min, max+step, step):

		fw_filename = sampleid + "_" + str(cnt) + "_FW.fastq"
		fwout = os.path.join(outdir, fw_filename)
		rv_filename = sampleid + "_" + str(cnt) + "_RV.fastq"
		rvout = os.path.join(outdir, rv_filename)

		seed = random.randint(100, 10000)

		if i > max:
			i = max

		cmd = "seqtk sample -s" + str(seed) + " " + fw_readfile + " " + str(i) + " > " + fwout
		cmd = cmd + " && "  "seqtk sample -s" + str(seed) + " " + rv_readfile + " " + str(i) + " > " + rvout
		data[cnt] = {'readcnt':i , 'cmd':cmd}
		cnt = cnt + 1
		#logger.debug("Cmd inside loop:"+ cmd)
		#fh.write(cmd + "\n")
	
	#fh.close()

	logger.info("Total cmds: {}".format(len(data)))	
	return data;

##########################################################################################
#
#
##########################################################################################

def srun_script(scriptfile):
        """This function will submit jobs to cluster"""
        logger.info("Submitting jobs to cluster")
        logger.debug("Inside function:" + scriptfile)
        cmd = "/scratch/gdlab/manish/scripts/ngs/submit_job_array_slurm.sh " + scriptfile + " 1 | sbatch --mem=30000M";
        proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
        output = proc.communicate()[0].strip()
        #output = subprocess.check_output(cmd, shell=True)
        logger.info(output)
        matchObj = re.search('^Submitted batch job (\d+)', output)
        if matchObj:
                jobId = matchObj.group(1)
                logger.debug("jobId = " + jobId)
        else:
                print "No match!"
                logger.exception("Error occured: No match!!. Could not get the job Id")
                sys.exit(1)

        status = wait_for_completion(jobId)
        return;

##########################################################################################

def wait_for_completion(jobId):
        """This function will halt the program untill all the jobs in an array get completed"""
        logger.info("Waiting for all jobs to finish")
        logger.debug("JobId inside function = " + jobId)
        flag = 1
        if jobId:
                while(flag):
                        time.sleep(5)
                        cmd = "squeue -j " + jobId + " -h -o %t"
                        #print "Command= ", cmd
                        job_state = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                        #logger.debug("job_state= {}".format(job_state))

                        if job_state =="R":
                                logger.debug("Job is running")
                        elif job_state == "PD":
                                logger.debug("Job is pending")
                        elif job_state == "CA" or job_state == "F" or job_state == "NF":
                                logger.exception("Job id :" + jobId + "failed. Exiting!!")
                                sys.exit(1)
                        elif job_state == "CD" or job_state == "":
                                logger.info("Job finished. Checking the status of all tasks!")
                                cmd = "sacct -P -n -j " + jobId + " --format jobid,jobname,elapsed,ReqMem,alloccpus,state,exitcode"
                                job_detail = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                                if job_detail:
                                        #logger.debug(job_detail)
                                        for line in job_detail.splitlines():
                                                logger.debug(line)
                                                status = line.split('|')
                                                logger.debug("job state " + status[5])
                                                if status[5] != "COMPLETED" or status[6]!="0:0":
                                                        logger.exception("Error Occurred...Exiting!!")
                                                        sys.exit(1)
                                else:
                                        logger.exception("Error occurred....Exiting!!")
                                        sys.exit(1)
                                flag = 0
        else:
                logger.exception("jobId is empty! Exiting!!")
                sys.exit(1)

        return;

#########################################################################################

if __name__ == "__main__":
	main(sys.argv[1:])
