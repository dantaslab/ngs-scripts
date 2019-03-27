#!/usr/bin/env python

"""
File Name       : parse_card_MB.py
Author          : Manish Boolchandani
Created On      : 2017-05-10
Last Modified   : 2019-03-22
Description     : A program to parse CARD ontology file (card.json and aro.obo) to get metadata
                  for the given list of ARO ids. The output file will contain the following 
                  information for each ARO id: 
                    1) AMR gene family: Gene family the resistance gene belongs to
                    2) AR Category: Higher level gropuping of AMR gene family
                    3) Resistace Mechanism: Mechanism by which a gene confers resistance
                    4) Drug Class: Class of drugs to which a gene confers resistance
                    5) Antibiotics: name of specific drugs
                    6) Efflux Component: whether a gene is part of an efflux pump
                    7) Efflux Regulator: whether a gene is a part of regulatory system for efflux proteins
                    8) Part_of_Efflux: lists the name of efflux pump the gene belongs to
                    9) Regulates: lists the name of efflux complex, the gene regulates 
                    10) NCBI Taxon ID:
                    11) NCBI Taxon Name:
                    12) Description:Gene description

Dependencies    : 
Usage           : parse_card_MB.py --infile aroids.txt --card_json card.json --aro_obo aro.obo --outfile card_mappingfile.txt 

CHANGE LOG      :

TODO            :
"""

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
logger = logging.getLogger(__name__)

def main(argv):

    # Parse arguments from the command line
    args = parse_arguments(argv)

    # Check the format of the required input files
    check_requirements(args)

    # Configure the logger
    log_file = "zz_parse_card_ontology.log"
    logging.basicConfig(filename=log_file,format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', 
                        level=getattr(logging, args.log_level), 
                        filemode='w', 
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    card_data = defaultdict(dict)
    final_output = defaultdict(dict)
    aroid_list = []

    if not args.out_file:
        out_file = os.path.splitext(os.path.basename(args.input_file))[0] + '_out.txt'
    else:
        out_file = args.out_file

    logger.info("Started pipeline to parse CARD ")
    logger.info("Input file = {}".format(args.input_file))
    logger.info("CARD json file = {}".format(args.card_file))
    logger.info("CARD ontolgy file = {}".format(args.ontology_file))

    # Read input file containing ARO ids
    aroid_list = read_inputfile(args.input_file, 1)
    logger.info("Total num of ARO ids in the input file = {}".format(str(len(aroid_list))))

    # Read CARD json file to get metadata for each ARO id
    card_data = read_card_json(args.card_file, aroid_list)

    if len(card_data):
        # Read ontology file in obo format to map relationship 
        # between efflux components and regulatory elements
        final_output = read_aro_obo(args.ontology_file, card_data)

        # Finally write the output to a file
        if len(final_output):
            write_tofile(out_file, final_output)
    else:
        print("Nothing to write..Exiting!")

    logger.info("Pipeline ended!!")


def parse_arguments(argv):
    parser = argparse.ArgumentParser( prog = 'parse_card_MB.py',
                description = 'A program to parse CARD database file in json format')

    parser.add_argument( '--version', action = 'version', 
                version = "%(prog)s v" + VERSION)
    
    parser.add_argument( '-i', '--infile', dest = "input_file", required=True,
                help = "Enter input file name that contains list of ARO ids")

    parser.add_argument( '--card_json', dest = "card_file", required = True,
                help = "path to card json file (card.json)")

    parser.add_argument( '--aro_obo', dest = 'ontology_file', required = True,
                help = "path to card ontology file (aro.abo)")

    parser.add_argument( '-o', '--outfile', dest = "out_file",
                help = "name of the output dir")

    parser.add_argument( '--log-level', dest = 'log_level', default = 'DEBUG',
                choices = ['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                help = 'level of messages to display in log file [DEFAULT: INFO]\n\n'
    )

    return parser.parse_args()


def check_requirements(args):
    input_files = [args.input_file, args.card_file, args.ontology_file]

    for infile in input_files:
        filepath = os.path.abspath(infile)
        try:
            s = os.stat(filepath)
            if s.st_size == 0:
                #logger.exception("Error Occurred! The file is either empty or does not exist = {}".format(filepath))
                sys.exit("CRITICAL ERROR: Can't read the input file selected = {}".format(filepath))
            else:
                return True
        except OSError as e:
            print("CRITICAL ERROR:The file does not exists: " + str(e))
            sys.exit(1)


def read_inputfile(filename, col):
    logger.info("Read Mapping file: " + filename)
    id_list = []
    col_index = int(col) - 1

    if col_index < 0:
        logger.error('Invalid column number is specified for read the input file')
        sys.exit(1)

    with open(filename) as fh:
        for line in fh:
            columns = line.rstrip('\n').split('\t')
            if len(columns) > int(col_index):
                id_list.append(columns[col_index])
            else:
                logger.error('Invalid column number is specified for read the input file')
                sys.exit(1)

    return id_list


def read_card_json(card_file, aroid_list):
    logger.info("Step 1: Read card.json file")    

    output_data = defaultdict(dict)

    # load card.json file to dictionary
    with open(card_file, 'r') as fh:
        card_data = json.load(fh)
  
    #print(next(iter(data.values()))) 

    for id in aroid_list:
        output_data[str(id)] = {}
        sample_id = id.rstrip('\n').split(':')[1]
        logger.debug("Sample id => {}".format(sample_id))

        for k,val in card_data.items():
            #logger.debug("{} => {}".format(str(k), str(val)))
            tmp_data = defaultdict(dict)

            if "ARO_accession" in val:
                aro_id = card_data[k]['ARO_accession']
                #logger.debug("aro id => {}".format(aro_id))

                if str(aro_id) == str(sample_id):
                    aro_name = card_data.get(k).get('ARO_name')
                    aro_desc = card_data.get(k).get('ARO_description')
                    output_data[str(id)]['aro_name'] = str(aro_name)
                    output_data[str(id)]['aro_desc'] = str(aro_desc.encode('utf-8'))
                    logger.info("{} => {}".format(str(id), str(aro_name)))

                    if val.get('model_sequences', {}).get('sequence'):
                        ncbi_model = val.get('model_sequences', {}).get('sequence')
                        for model_id, model_val in ncbi_model.items():
                            if "NCBI_taxonomy" in model_val:
                                taxon_id = model_val['NCBI_taxonomy'].get('NCBI_taxonomy_id','NA')
                                taxon_name = model_val['NCBI_taxonomy'].get('NCBI_taxonomy_name','NA')
                                tmp_data['taxon_info'].update({str(taxon_id):str(taxon_name)})
                                #tmp_data['taxon_id'].append(str(taxon_id))
                                #tmp_data['taxon_name'].append(str(taxon_name))


                    if "ARO_category" in val:
                        aro_categories = card_data.get(k).get('ARO_category','')
                        #logger.debug("{} => {}".format(str(id),str(aro_categories)))

                        if aro_categories:
                            for catid, catval in aro_categories.items():
                                class_name = catval.get('category_aro_class_name')

                                if class_name.lower() == "resistance mechanism":
                                    res_acc = "ARO:" + catval.get('category_aro_accession')
                                    res_name = catval.get('category_aro_name')
                                    tmp_data['res_mechanism'].update({str(res_acc):str(res_name)})
                                    #tmp_data['res_mechanism'].append(str(res_name))

                                elif class_name.lower() == "amr gene family":
                                    argfam_acc = "ARO:" + catval.get('category_aro_accession')
                                    argfam_name = catval.get('category_aro_name')
                                    tmp_data['arg_family'].update({str(argfam_acc):str(argfam_name)})
                                    #tmp_data['arg_family'].append(str(argfam_name))

                                elif class_name.lower() == "antibiotic":
                                    abx_acc = "ARO:" + catval.get('category_aro_accession')
                                    abx_name = catval.get('category_aro_name')
                                    tmp_data['antibiotics'].update({str(abx_acc):str(abx_name)})
                                    #tmp_data['antibiotics'].append(str(abx_name))

                                elif class_name.lower() == "drug class":
                                    drug_acc = "ARO:" + catval.get('category_aro_accession')
                                    drug_name = catval.get('category_aro_name')
                                    tmp_data['drug_class'].update({str(drug_acc):str(drug_name)})
                                    #tmp_data['drug_class'].append(str(drug_name))
                                
                                elif class_name.lower() == "efflux component":
                                    efflux_acc = "ARO:" + catval.get('category_aro_accession')
                                    efflux_name = catval.get('category_aro_name')
                                    tmp_data['efflux_component'].update({str(efflux_acc):"Yes"})
                                    #tmp_data['efflux_componenet'] = "Yes"

                                elif class_name.lower() == "efflux regulator":
                                    reg_acc = "ARO:" + catval.get('category_aro_accession')
                                    reg_name = catval.get('category_aro_name')
                                    tmp_data['regulator'].update({str(reg_acc):"Yes"})
                                    #tmp_data['regulator'] = "Yes"                                    
                                         
                                else:
                                    logger.info("class_name => {} not included".format(str(class_name))) 


                            output_data[str(id)].update(tmp_data)
                            logger.debug("output stream = {}".format(str(output_data[str(id)])))
    
    return output_data


def read_aro_obo(ontology_file, card_data):
    """
    Read ontology file and detect efflux pump
    """
    logger.info("Step 2: Read ontology file") 

    efflux_data = defaultdict(dict)
    eflx_ids = []
    reg_ids = []
    uniq_argfam = []
    amr_families = defaultdict(list)

    for aroid, data in card_data.iteritems():
        arg_fams = data.get('arg_family', {})

        eflx_aro = data.get('efflux_component', '')
        reg_aro = data.get('regulator','')
    
        for famid, fam in arg_fams.items():
            amr_families[str(aroid)].append(str(famid))

        if eflx_aro and eflx_aro != 'NA':
            eflx_ids.append(aroid)

        if reg_aro and reg_aro != 'NA':
            reg_ids.append(aroid)

    #print("card_data = {}\n".format(str(card_data)))
    #print("amr_families = {}\n".format(str(amr_families)))

    #path = getpath(amr_families, 'ARO:3000237')
    #print("path = {}".format(str(path)))
    
    logger.info("Num of unique AMR gene family = {}".format(str(len(amr_families))))
    logger.info("Num of efflux genes = {}".format(str(len(eflx_ids))))
    logger.info("Num of regulator genes = {}".format(str(len(reg_ids))))
    parent_data = defaultdict(dict) 
    aro_onto = defaultdict(dict)

    with open(ontology_file, 'r') as f:
        for i, group in enumerate(get_groups(f, "[Term]\n"), start=1):

            aro_id = group.get('id',[])
            #print("aroid = {}".format(str(aro_id)))
            aro_name = group.get('name',[])
            eflx_comp = []
            data = defaultdict(list)
            parent_terms = group.get('is_a', [])


            if len(aro_id)== 1 and aro_id[0].startswith('ARO:') and parent_terms:
                aro_onto[str(aro_id[0])] = {'aro_name':'NA', 'parent_terms':{}, 'ontology':'NA'}
                aro_onto[aro_id[0]]['aro_name'] = str(aro_name[0])
                term_list = defaultdict(list)
                for item in parent_terms:
                    terms = re.findall('^(ARO:\d+)\s!\s(.*)$', item)
                    if terms:
                        term_list['parent_terms'].append(terms[0][0])

                aro_onto[str(aro_id[0])].update(term_list)

            #if len(aro_id) == 1 and aro_id[0] in [x for v in amr_families.values() for x in v]:
                #parent_terms = group.get('is_a', [])
                #for k, v in amr_families.items():
                #    for idx, val in enumerate(v):
                #        if val == aro_id[0]:
                #            print(k + "=>" + str(idx) + "=>" + val)
                
            for k, v in amr_families.items():
                tmp_data = defaultdict(dict)
                for idx, val in enumerate(v):
                    if len(aro_id) == 1 and val == aro_id[0]:
                        #parent_terms = group.get('is_a', [])
                        for item in parent_terms:
                            tmp = re.findall('^(ARO:\d+)\s!\s(.*)$', item)
                            tmp_data[val].update({tmp[0][0]:tmp[0][1]})
                            #print(k + "=>" + str(tmp_data))

                parent_data[k].update(tmp_data)
                #logger.debug("{} => {}\n".format(str(aro_id[0]), str(data)))

            if len(aro_id) == 1 and aro_id[0] in eflx_ids:
                #logger.debug("Record #{} => {}\n".format(str(aro_id), str(group)))
                aro_name = group.get('name',[])
                rel = group.get('relationship',[])

                for line in rel:
                    eflx_comp = re.findall('^part_of\s(ARO:\d+)\s!\s(.*)$', line)
                    reg_info = re.findall('^regulates\s(ARO:\d+)\s!\s(.*)$', line)

                    if eflx_comp:
                        eflx_comp_id = eflx_comp[0][0]
                        eflx_comp_name = eflx_comp[0][1]
                        data['efflux_pump'].append(eflx_comp_name)

                    if reg_info:
                        reg_id = reg_info[0][0]
                        reg_name = reg_info[0][1]
                        data['regulates'].append(str(reg_name))
                
                logger.debug("{} => {}".format(aro_name, str(data)))
                card_data[str(aro_id[0])].update(data)

    #print("parent_data = {}\n".format(str(parent_data)))
    print("Num of aro_data = {}\n".format(str(len(aro_onto))))

    for k, v in aro_onto.items():
        ontology = []
        parents = v.get('parent_terms',[])
        ontology.append("0-" + k)
        for idx, item in enumerate(parents, 1):
            ontology.append("{}-{}".format(str(idx), str(item)))
            count = 2
            get_aro_parents(item, aro_onto, ontology, count)

        ontology = list(filter(None, ontology))
    
        aro_onto[k].update({'ontology' : '|'.join(ontology)})


    #print("aro_onto = {}\n".format(str(aro_onto)))
    #keys = [k for k,v in aro_onto.items() if 'ARO:3000014' in v]
    #print('key = {}'.format(str(len(keys))))

    amr_category_data = defaultdict(dict)
    final_data = defaultdict(dict)

    for k in card_data.keys():
        cat = parent_data.get(k, {})
        ont = aro_onto.get(k, {})
        print("ont = {}".format(str(ont)))
        print("name = {}".format(ont.get('aro_name','NA')))
        amr_category_data[str(k)] = {'arg_category': {}}
        if cat:
            for x, v in cat.items():
                amr_category_data[str(k)]['arg_category'].update(v)
                #print("{} => {}".format(str(k), str(v)))
            card_data[k].update(amr_category_data[k])

        if ont:
            tmp = ont.get('ontology','')
            #print("{}={}".format(str(k),str(tmp)))
            ont_list = []
            for x in sorted(tmp.split('|')):
                y = str(x.split('-')[1])
                z = int(x.split('-')[0])
                if y and z > 1: 
                    print("{} => {}".format(str(y), aro_onto[str(y)].get('aro_name','NA')))
                    ont_list.append(aro_onto[str(y)].get('aro_name','NA'))

            card_data[k].update({'arg_ontology':ont_list})

    #print("arg_category = {}\n".format(str(amr_category_data)))
    #print("card_data = {}\n".format(str(card_data)))

    return card_data

def get_aro_parents(aro_parent, data, ontology, count):
    #print(str(count) + " > " + str(ontology))
    for k,v in data.items():
        if aro_parent == k:
            parents = v.get('parent_terms',[])
            for x in parents:
                if not any(x in i for i in ontology):
                    ontology.append(str(count)+"-"+x)
            count+=1
            for i in parents:
                if i.startswith('ARO:') and i != 'ARO:3000000':
                    get_aro_parents(i, data, ontology, count)



def write_tofile(outfile, final_output):
    """
    Write final output to tab delimited file
    """
    logger.info("Step 3: Write output file = {}".format(str(outfile)))

    if len(final_output):
        with open(outfile, 'w') as f1:
            f1.write("ARO_Id\tARO_Name\tAMR_Gene_Family\tAMR_Category\tAMR_Ontology\tResistance_Mechanism"\
                    "\tDrug_Class\tAntibiotics\tEfflux_Component\tEfflux_Regulator"\
                    "\tPart_of_Efflux_Pump\tRegulates"\
                    "\tNCBI_Taxon_Id\tNCBI_Taxon_Name\tDescription\n")

            for id, val in final_output.iteritems():
                sample_id = str(id)
                tmp = {}

                for k, v in val.iteritems():
                    if isinstance(v, list):
                        tmp[str(k)] = '|'.join(map(str, v))
                    elif isinstance(v, dict):
                        tmp[str(k)] = '|'.join(map(str, v.values()))
                    else:
                        tmp[str(k)] = str(v)
    
                f1.write("{id}\t{aro_name}\t{arg_fam}\t{arg_cat}\t{arg_ont}\t{res_mec}\t{abx_class}\t{abx_name}"\
                        "\t{efflux}\t{regulator}\t{pump}\t{regulates}\t{taxid}\t{taxname}\t{desc}\n"\
                        .format(
                                id = sample_id,
                                aro_name = tmp.get('aro_name','NA'),
                                arg_fam = tmp.get('arg_family', 'NA'),
                                arg_cat = tmp.get('arg_category','NA'),
                                arg_ont = tmp.get('arg_ontology','NA'),
                                res_mec = tmp.get('res_mechanism', 'NA'),
                                abx_class = tmp.get('drug_class', 'NA'),
                                abx_name = tmp.get('antibiotics', 'NA'),
                                efflux = tmp.get('efflux_component','NA'),
                                regulator = tmp.get('regulator','NA'),
                                pump = tmp.get('efflux_pump', 'NA'),
                                regulates = tmp.get('regulates','NA'),
                                taxid = tmp.get('taxon_id','NA'),
                                taxname = tmp.get('taxon_name','NA'),
                                desc = tmp.get('aro_desc', 'NA')
                         ))



def file_exists(filepath):
    """ 
    Check if the given file exists and non empty
    """
    try:
        s = os.stat(filepath)
        if s.st_size > 0:
            return True
        else:
            logger.exception("Error Occurred! The file {} is empty".format(ontology_file))
            return False
    except OSError as e:
        logger.exception("Error Occurred: " + str(e))
        sys.exit(1)


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
    data = defaultdict(list)
    for line in seq:
        # Here the `startswith()` logic can be replaced with other
        # condition(s) depending on the requirement.
        if line.startswith(group_by):
            if data:
                yield data
                data = defaultdict(list)

        if ': ' in line:
            tmp = line.rstrip('\n').split(': ')
            key = tmp[0].strip()
            val = tmp[1].strip()
            data[key].append(str(val))

    if data:
        yield data


if  __name__ == "__main__":
    main(sys.argv[1:])
