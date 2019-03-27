#!/usr/bin/env python

'''
File Name	: parse_card_json.py
Author		: Manish Boolchandani
Created On	: 2017-05-10
Description	: A program to parse CARR json file

'''
import sys

import os
import operator
import time
import re
import json

import logging
import logging.handlers
import argparse
import subprocess
from collections import defaultdict
from itertools import groupby

VERSION = "1.0.0"

#Logger configuration and setup
LEVEL = logging.DEBUG
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y%I:%M:%S %p')
logger = logging.getLogger(__file__)
logger.setLevel(LEVEL)
handler = logging.handlers.RotatingFileHandler('zz_parseJSON.log', maxBytes=2000000, backupCount=20)
handler.setLevel(LEVEL)
handler.setFormatter(formatter)
logger.addHandler(handler)


def main(argv):
	logger.info("A program to parse CARD json file")
	
	args = parse_arguments(argv)
	
	id_list = read_inputfile(args.input_file, 1)
	
	if not args.out_file:
		out_file = os.path.splitext(os.path.basename(args.input_file))[0] + '_out.txt'
	else:
		out_file = args.out_file

        #logger.info("IDs: {}".format(str(id_list)))
	#data = read_card_json(args.card_file)
	
	logger.info("input file = {}".format(args.input_file))
	logger.info("CARD database file = {}".format(args.card_file))
	logger.info("Total num of ARO ids = {}".format(str(len(id_list))))
	
	output_data = defaultdict(dict)
	
	with open('input.txt') as f:
    		for i, group in enumerate(get_groups(f, "[Term]\n"), start=1):
        		print ("Record #{}".format(i))
			aro_id = None
			aro_name = None
			parent_term=[]
			relationship = []

			for line in group:
				if line.startswith('id'):
					aro_id = line.rstrip('\n').split('id: ')[1]
					print "id => {}".format(str(aro_id))
				
				if str(aro_id):		
					
					output_data[str(aro_id)] = {'aro_name':'NA', 'parent_terms':'NA', 'relationship':'NA'}
					
					if line.startswith('name:'):
						aro_name = line.rstrip('\n').split('name: ')[1]
						output_data[str(aro_id)]['aro_name'] = str(aro_name)
						print "name => {}".format(str(aro_name))
					if line.startswith('is_a:'):
						parent_term.append(line.rstrip('\n').split('is_a: ')[1])
						output_data[str(aro_id)]['parent_terms'] = '|'.join(parent_term)

						print "parent => {}".format(str(parent_term))
				
					if line.startswith('relationship:'):
						relationship.append(line.rstrip('\n').split('relationship: ')[1])
						output_data[str(aro_id)]['relationship'] = '|'.join(relationship)

						print "relationship => {}".format(str(relationship))


        		#print ("".join(group))
	
	#with open('input.txt') as f:
    	#	for k, group in groupby(f, key=make_grouper("[Term\n]")):
        #		# parse >header description
        #		header, description = next(group)[1:].split(maxsplit=1)
        #		for line in group:
        #    			# handle the contents of the section line by line
	#			print "line => {}".format(str(line))	
	
	#for id in id_list:
	#	output_data[str(id)] = {'aro_name':'NA', 'category_aro_id':'NA', 'category_aro_name':'NA'}
	#	sample_id = id.split(':')[1]
	#	#item_vals = getpath(data, str_id)
	#	for k,val in data.items():
	#		#logger.info("{} => {}".format(str(k), str(val)))
	#		if "ARO_accession" in val: 
	#			aro_id = data[k]['ARO_accession']
	#			if str(aro_id) == str(sample_id):
	#				aro_name = data.get(k).get('ARO_name')
	#				output_data[str(id)]['aro_name'] = str(aro_name)
	#				logger.info("{} => {}".format(str(id),str(aro_name)))
	#				if "ARO_category" in val:
	#					aro_categories = data.get(k).get('ARO_category','')
	#					#logger.info("{} => {}".format(str(id),str(aro_categories)))
	#					if aro_categories:
	#						all_cat_aro = []
	#						all_cat_name = []
	#						for catid, catval in aro_categories.items():
	#							cat_aro_acc = catval.get('category_aro_accession')
	#							cat_aro_name = catval.get('category_aro_name')
	#							all_cat_aro.append(str(cat_aro_acc))
	#							all_cat_name.append(str(cat_aro_name))
	#							#logger.info("{} => {}".format(cat_aro_acc,cat_aro_name))
	#							
	#						output_data[str(id)]['category_aro_id'] = '|'.join(all_cat_aro)
	#						output_data[str(id)]['category_aro_name'] = '|'.join(all_cat_name) 
	#	
	#
	logger.info("output stream = {}".format(str(output_data)))
	#write_tofile(out_file, output_data)


def parse_arguments(argv):

        parser = argparse.ArgumentParser(
                prog = 'parse_card_json.py',
                description = 'A program to parse CARD database file in json format')
        parser.add_argument(
                '--version',
                action = 'version',
                version = "%(prog)s v" + VERSION)
        parser.add_argument(
                '-i', '--infile',
                dest = "input_file",
                required=True,
                help = "Enter input file name that contains list of ARO ids")
        parser.add_argument(
                '-d', '--card',
                dest = "card_file",
                required = True,
                help = "path to sequence directory")
        parser.add_argument(
                '-o', '--outfile',
                dest = "out_file",
                help = "name of the output dir")
        
        return parser.parse_args()


def read_inputfile(filename, col):
	logger.info("Reading Mapping file: " + filename)
	id_list = []
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
                               id_list.append(columns[col_index])
                        else:
                                logger.error('Invalid column number is specified for read the input file')
                                sys.exit(1)

        return id_list

def read_card_json(card_file):
	
	try:
        	s = os.stat(card_file)
                if s.st_size == 0:
                        logger.exception("Error Occurred! The file {} is empty".format(card_file))
                        sys.exit(1)
        except OSError as e:
                logger.exception("Error Occurred: " + str(e))
                sys.exit(1)
	
	with open(card_file, 'r') as fh:
		data = json.load(fh)

	return data

def write_tofile(outfile,output_data):
	
	with open(outfile, 'w') as f1:
		for id, val in output_data.iteritems():
			sample_id = str(id)
			aro_name = val['aro_name']
			aro_category_id = val['category_aro_id']
			aro_category_name = val['category_aro_name']
			f1.write("{}\t{}\t{}\t{}\n".format(sample_id, aro_name, aro_category_id, aro_category_name)) 


def find_key(d, value):
    for k,v in d.items():
        if isinstance(v, dict):
            p = find_key(v, value)
            if p:
                return [k] + p
        elif v == value:
            return [k]

def getpath(nested_dict, value, prepath=()):
    for k, v in nested_dict.items():
        path = prepath + (k,)
        if v == value: # found value
            return path
        elif hasattr(v, 'items'): # v is a dict
            p = getpath(v, value, path) # recursive call
            if p is not None:
                return p

def get_groups(seq, group_by):
    data = []
    for line in seq:
        # Here the `startswith()` logic can be replaced with other
        # condition(s) depending on the requirement.
        if line.startswith(group_by):
            if data:
                yield data
                data = []
        data.append(line)

    if data:
        yield data

#def make_grouper(group_by):
#    counter = 0
#    def key(line):
#        nonlocal counter
#        if line.startswith(group_by):
#            counter += 1
#        return counter
#    return key

if  __name__ == "__main__":
	main(sys.argv[1:])
