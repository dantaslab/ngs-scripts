#!/usr/bin/env python

"""
File Name    : combine_reads_MB.py
Author       : Manish Boolchandani, manish@wustl.edu
Created On   : 2016-03-28
Last Modified: 2016-03-28

Description  :  This program will combine read files of the same sample
                obtained from multiple sequencing runs


Dependencies:

Usage:  combine_reads_MB.py [-h] [-v] [-i INPUT_FILE] [-d SEQ_DIR]
                           [-o OUTPUT_DIR] [-m MIN_READCNT] [-x MAX_READCNT]
                           [-r TOTAL_READCNT] [-w]
CHANGE LOG:
TODO:

"""

# Python imports
import os, operator, sys, time, re, logging, argparse, subprocess, glob
import logging.handlers
from collections import defaultdict

# Version
VERSION = "1.0"
# Logger configuration and setup
LEVEL = logging.DEBUG
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger(__file__)
logger.setLevel(LEVEL)
handler = logging.handlers.RotatingFileHandler('zz_mergeReads.log', maxBytes=2000000, backupCount=20)
handler.setLevel(LEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)

def main(argv):
    
    logger.info("STARTED PIPELINE!!")
    parser = argparse.ArgumentParser(prog = 'combine_reads_MB.py', description = "A program to combine reads of a sample from multiple sequencing runs")
    parser.add_argument('-v','--version', action = 'version', version = "%(prog)s v" + str(VERSION))
    parser.add_argument('-i', '--ifile', dest='input_file', help="input file containing sample ids")
    parser.add_argument('-d', '--seqdir', dest='seq_dir', help="path to sequence directory")
    parser.add_argument('-o','--outdir', dest='output_dir', help="Name of an output directory")
    parser.add_argument('-m','--min_readcnt', dest='min_readcnt', help="Enter minimum reads to consider while merging files")
    parser.add_argument('-x', '--max_readcnt', dest='max_readcnt', help="Maximun number of reads per sample")
    parser.add_argument('-r','--total_readcnt', dest='total_readcnt', help="Enter total read count you want to keep in the output directory")
    parser.add_argument('-w','--overwrite', action='store_true', default='False', help="overwrite, if the output files exists")
    args = parser.parse_args()

    # check to see that input file is provided
    if not args.input_file:
        logger.exception("You must provide input file containing sample ids.")
        parser.exit(status=1, message="You must provide input file containing sample ids\n\n")
    else:
        input_file = args.input_file

    if not args.output_dir:
        outputdirname = os.path.splitext(os.path.basename(input_file))[0] + "_out"
        output_dir = os.path.join(os.getcwd(), outputdirname)
        logger.warning("Output directory is not provided. Writing the result to {} " . format(output_dir))
    else:
        output_dir = os.path.realpath(args.output_dir)

    if not args.seq_dir:
        logger.exception("You must provide path to directory that contains sequences or metaphlan output")
        parser.exit(status=1, message="You must provide path to sequence directory\n")
    else:
        seq_dir = os.path.realpath(args.seq_dir)

    if not args.min_readcnt:
        min_readcnt = 100000
    else:
        min_readcnt = int(args.min_readcnt)

    if not args.max_readcnt:
        max_readcnt = 5000000    
    else:
        max_readcnt = int(args.max_readcnt)

    if not args.total_readcnt:
        readcnt_tomerge = 4000000
    else:
        readcnt_tomerge = int(args.total_readcnt)

    readcnt_millions = str(readcnt_tomerge/float(1000000))
    #print "readcnt in millions {}".format(str(readcnt_millions))
    file_suffix = str(readcnt_millions + 'M.txt')

    try:
        if not os.path.isdir(output_dir): subprocess.check_call(['mkdir', output_dir])
        datadict = read_inputfile(input_file)
        sample_stat = {}
        #logger.info("data = {}".format(str(datadict)))
        fh1 = open('r01_combined.txt', 'w')
        fh2 = open('r02_sampleids_gt' + str(file_suffix), 'w')
        fh3 = open('r03_sampleids_lt' + str(file_suffix), 'w')
        fh1.write('#SampleLabNumber\tSampleSeqID\tReadCnt\n')       
        fh2.write('#SampleLabNumber\tSampleSeqID\tReadCnt\n')
        fh3.write('#SampleLabNumber\tSampleSeqID\tReadCnt\n')       
        
        for id, data in datadict.iteritems():
        
            total_cnt = 0;
            seqids = ""
            #logger.info("id = " + str(id))

            for seqid, readcnt in data.iteritems():
                if int(readcnt) > int(min_readcnt):
                    total_cnt += int(readcnt)
                    seqids += str(seqid) + ";"  
                    #total_cnt = sum(int(v) for v in data.values() if int(v) > 100000)
                    #seqids = ";".join([str(x) for x in data.keys() if data[x] > 100000])

            if len(seqids):             
                seqids = seqids.rstrip(';')
                fh1.write(id + '\t' + seqids + '\t' + str(total_cnt) + '\n')
        
                if total_cnt >= int(readcnt_tomerge):
                    seqid_list = seqids.split(';')
                    participant_id = seqid_list[0].split('_')[0]
                    sample_id = str(participant_id + "_" + seqid_list[0].split('_')[1])
                    fw_cmd = ''
                    rv_cmd = '' 
                    
                    if len(seqid_list) > 1:
                        fw_file = ''
                        rv_file = ''

                        for seqid in seqid_list:        
                            fw_file += os.path.join(seq_dir, str(seqid + '_FW_CLEAN.fastq')) + ' '
                            rv_file += os.path.join(seq_dir, str(seqid + '_RV_CLEAN.fastq')) + ' '
                        
                        fw_outfile = os.path.join(output_dir, str(sample_id + '_COMB_FW_CLEAN.fastq'))
                        rv_outfile = os.path.join(output_dir, str(sample_id + '_COMB_RV_CLEAN.fastq'))
                        fh2.write(str(id) + '\t' + str(sample_id +'_COMB') + '\t' + str(total_cnt) + '\n')
                    else:
                        fw_file = ''
                        rv_file = ''
                        fw_file = os.path.join(seq_dir, str(seqid_list[0] + '_FW_CLEAN.fastq'))
                        rv_file = os.path.join(seq_dir, str(seqid_list[0] + '_RV_CLEAN.fastq'))
                        fw_outfile = os.path.join(output_dir, str(seqid_list[0] + '_FW_CLEAN.fastq'))
                        rv_outfile = os.path.join(output_dir, str(seqid_list[0] + '_RV_CLEAN.fastq'))
                        fh2.write(str(id) + '\t' + str(seqid_list[0]) + '\t' + str(total_cnt) + '\n')
                    

                    if not os.path.isfile(fw_outfile):
                        fw_cmd = 'cat {} > {}'.format(fw_file, fw_outfile)
                        rv_cmd = 'cat {} > {}'.format(rv_file, rv_outfile)
                        logger.debug("fw_cmd = {} and rv_cmd = {}".format(fw_cmd, rv_cmd))
                        if fw_cmd and rv_cmd:
                            proc1 = subprocess.Popen(fw_cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
                            output1 = proc1.communicate()[0].strip()
                            proc2 = subprocess.Popen(rv_cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
                            output2 = proc2.communicate()[0].strip()
                            #pass
                else:
                    fh3.write(id + '\t' + seqids + '\t' + str(total_cnt) + '\n')


                #if int(total_cnt) > int(max_readcnt):
                #     logger.info("id => {} has more than {} reads".format(str(id), str(max_readcnt)))           
        
        fh1.close()
        fh2.close()
        fh3.close() 
        

    except:
        logger.exception("Error Occurred")
        sys.exit(1)

    logger.info("End of pipeline!!")


def read_inputfile(filename):
    logger.info("Reading input file: " + filename)
    sample_info = defaultdict(dict)
    try:
        s = os.stat(filename)
        if s.st_size == 0:
            logger.exception("Error Occurred! The file {} is empty".format(filename))
            sys.exit(1)
    except OSError as e:
        logger.exception("Error Occurred: " + str(e))
        sys.exit(2)

    fh = open(filename, "r")
    data = [line.rstrip('\n').split("\t") for line in fh if not line.startswith('#')]
    fh.close()
    #logger.debug("data = {}".format(str(data)))
    for item in data:
        sample_info[item[0]][item[1]] = int(item[2].replace(',',''))
    
    logger.info("IDs: {}".format(str(len(sample_info))))

    return sample_info

####################################################################

if __name__ == '__main__':
    main(sys.argv[1:])
