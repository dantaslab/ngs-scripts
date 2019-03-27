#!/usr/bin/env python

"""
File Name    :  preprocess_config_MB.py
Author       :  Manish Boolchandani, manish@wustl.edu
Created On   :  2017-04-04
Last Modified:  2017-04-25

Description :   This is a program to parse config file and set values 
                that will be used in the subsequent programs. The config
                file "preprocess_seqreads.cfg" will store parameters and their 
                values for various programs like trimmomatic, deconseq etc.

"""

# Python imports
import sys
import os

import argparse
import logging
import ConfigParser

# Logger configuration and setup
logger = logging.getLogger(__name__)

# Config file path
cnfpath = "/scratch/gdlab/manish/travelersDiarrhea/scripts/preprocess_config.cfg"

# log options
log_level_choices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
log_level = log_level_choices[1]

class PreprocessConfig:
	def __init__(self,trimmomatic, deconseq, general):
		self.trimmomatic = trimmomatic
		self.deconseq = deconseq
		self.general = general

def read_configfile(cnfpath):
	# get the settings from config file
	config_items = parse_configfile(cnfpath)
	
	# get general parameter values
	f_suffix = get_item(config_items,"general", "file_suffix_fw", "string")
	r_suffix = get_item(config_items,"general", "file_suffix_rv", "string")
	general = {'file_suffix_fw':f_suffix, 
                   'file_suffix_rv':r_suffix}
	
	# quality trimming using trimmomatic
	illuminaclip = get_item(config_items, "quality_trimming", "illuminaclip", "string")
	adapt = get_item(config_items,"quality_trimming","adapt","string")
	leading = get_item(config_items, "quality_trimming", "leading", "int")
	trailing = get_item(config_items, "quality_trimming", "trailing", "int")
	sliding_window = get_item(config_items, "quality_trimming", "slidingwindow", "string")
	minlen = get_item(config_items, "quality_trimming", "minlen", "int")
	trim_out_suffix = get_item(config_items,"quality_trimming","file_suffix","suffix")
	
	trimmomatic = {'illuminaclip':illuminaclip,
			'adapt':adapt,
			'leading':leading,
			'trailing':trailing,
			'sliding_window':sliding_window,
			'minlen':minlen,
			'trim_out_suffix':trim_out_suffix	
			} 
	
	# deconseq options for host contamination removal
	#deconseq_path = get_item(config_items, "filtering", "deconseq_path", "string")
	#deconseq_path = "/home/manish/software/deconseq/0.4.3/deconseq.pl"
	deconseq_path = "/opt/apps/labs/gdlab/software/deconseq/0.4.3-chr38/deconseq.pl"
	deconseq_db = get_item(config_items,"filtering","deconseq_db","strings")
	deconseq_out_suffix = get_item(config_items,"filtering","file_suffix","string")
 
	deconseq = {'deconseq_path': deconseq_path, 'deconseq_db':deconseq_db, 'deconseq_out_suffix': deconseq_out_suffix}

	#Write to the log file the configuration used for the run
	loginf = []
	
	loginf.append("\n" + "-"*75)
	loginf.append("Config file = " + str(cnfpath))
	loginf.append("log level = " + log_level)
	
	loginf.append("-"*75)
	loginf.append("TRIMMING READS USING TRIMMOMMATIC")	
	loginf.append("")
	loginf.append("Adapter = " + str(adapt))
	loginf.append("Illuminaclip = " + str(illuminaclip))
	loginf.append("Leading = " + str(leading))
	loginf.append("Traling = " + str(trailing))
	loginf.append("Sliding Window =" + str(sliding_window))
	loginf.append("Min Length = " + str(minlen))
	
	loginf.append("-"*75)
	loginf.append("FILTERING READS USING DECONSEQ")
	loginf.append("")
	loginf.append("Deconseq path = " + str(deconseq_path))
	loginf.append("Deconseq database = " + str(deconseq_db))
	loginf.append("Deconseq output file suffix = " + str(deconseq_out_suffix))

	logger.info("\n\nRun Config Settings:" + "\n".join(loginf))
	loginf.append("-"*75)
	
	config_options = PreprocessConfig(trimmomatic, deconseq, general)

	return config_options

# Reading Config File
def parse_configfile(cnfpath):
        
	config = ConfigParser.ConfigParser()
        
	try:
                config.read(cnfpath)
        except OSError as e:
                logger.exception("Error: Unable to read the config file: " + str(e))

        config_list = {}

        for section_name in config.sections():
                #logger.info("Section:" + section_name)
                config_list[section_name] = {}
                section_list = config.items(section_name)
                #logger.info("Options:" + str(section_list))
                for name,value in section_list:
                       #logger.info("Section -> {}, Option -> {}, Value -> {}".format(section_name,name,value))
                       config_list[section_name][name] = value

        logger.debug("config list = " + str(config_list))
	return config_list



# Get the item from the dictionary of section/names from the user edit config file
def get_item(config_items, section, name, type=None):
    # try to obtain the value from the config dictionary
    try:
        value=config_items[section][name]
    except KeyError:
        sys.exit("CRITICAL ERROR: Unable to load config file" \
            " \nItem not found. \nItem should be in section (" + section + ") with name (" + name + ").")

    # if present, try to change the value type
    if type:
        try:
            if type == "string":
                value=str(value)
            elif type == "int":
                value=int(value)
            elif type == "float":
                value=float(value)
            elif type == "bool":
                if value in ["False","false","F","f"]:
                    value=False
                elif value in ["True","true","T","t"]:
                    value=True
                else:
                    raise ValueError
        except ValueError:
            sys.exit("CRITICAL ERROR: Unable to load value from config file" +
                     " \nItem found in section (" + section + ") with name (" + name + "). " +
                     "\nItem is not of type (" + type + ").")

    return value


