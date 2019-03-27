#!/usr/bin/env python

"""
File Name    : hmmscan-parser_MB.py
Author       : Manish Boolchandani, manish@wustl.edu
Created On   : 2018-05-01
Last Modified: 2018-05-01
Description  : Parses hmm output file (e.g Resfams-targets.txt) 
               and convert it to tabular format

Usage: hmmscan-parser_MB.py <hmm_outfile> 

CHANGELOG:
TODO:

"""

import os,sys
import collections

def main(args):
    hmm_file = args[0]
    hmm_dict = None
    f_name = os.path.splitext(os.path.basename(hmm_file))[0]
    
    if os.stat(os.path.abspath(hmm_file)).st_size > 0:
        hmm_dict = hmmscanParser(hmm_file)
        if hmm_dict:
            write_hmm_output(hmm_dict, f_name)
        else:
            print("Error while parsing the file!")
            sys.exit(1)
    else:
        print("Input hmm file is empty or doesnot exist!")
        sys.exit(1)

def hmmscanParser(hmm_file):
    hmm_dict = collections.defaultdict(dict)
    records = {}
    header = ["target_name", "target_accession", "query_name", "query_accession", "seq_eval",
                "seq_score", "seq_bias", "dom_eval", "dom_score","dom_bias",
                "exp_dom","num_reg", "num_clu", "num_overlap", "num_env",
                "num_dom", "num_rep", "num_inc", "description_of_target"]

    with open(hmm_file, 'r') as fhandle:
        for line in fhandle:
            rec = {}
            record = {}

            if line.startswith("#") or line is "":
                continue
            line = line.rstrip()
            itemlist = line.split(None, 18)
            for idx, item in enumerate(itemlist):
                rec[header[idx]] = item

            hmm_dict[itemlist[2]][itemlist[1]] = rec

    return hmm_dict

def write_hmm_output(output_dict, f_name):
        output_file = f_name + "_out.txt"
        if output_dict is not None:
            with open(output_file, 'w') as OUT:
                #header = ["target_name", "target_accession", "query_name", "query_accession", "seq_eval",
                #"seq_score", "seq_bias", "dom_eval", "dom_score","dom_bias",
                #"exp_dom","num_reg", "num_clu", "num_overlap", "num_env",
                #"num_dom", "num_rep", "num_inc", "description_of_target"]

                OUT.write("SampleId\tQueryName\tQueryAccession\tTargetName\tTargetAccession\tTargetDescription\tSeqEvalue\tSeqScore\n")
                for qid, records in output_dict.items():
                    for hitid, result in records.items():
                        OUT.write("{sample}\t{qname}\t{qacc}\t{sname}\t{sacc}\t{sdesc}\t{evalue}\t{score}\n".format(
                            sample = str(f_name),
                            qname = str(result["query_name"]),
                            qacc = str(result["query_accession"]),
                            sname = str(result["target_name"]),
                            sacc = str(result["target_accession"]),
                            sdesc = str(result["description_of_target"]),
                            evalue = str(result["seq_eval"]),
                            score = str(result["seq_score"]),
                        ))
        else:
            print("No result for {}".format(str(self.db_params['dbname'])))

if __name__=="__main__":
    main(sys.argv[1:])
