# ngs-scripts module : Set of scripts to process high-throughput sequencing data used in [Dantas Lab](http://www.dantaslab.org/)

#### __General utility scripts__
- __Preprocessing of sequencing data__
  - *preprocess_seqreads_MB.py* : Performs preprocessing of raw sequencing reads and gives the clean reads after quality trimming and removing host contaminants
  - *preprocess_config_MB.py*
  - *utilities_MB.py*
  - *check_job_completion_MB.py*
- __Other useful scripts__
  - *subsample_PEreads_MB.py* : sub-sampling of paired-end sequencing reads
  - *combine_reads_MB.py* : combine reads from multiple sequencing runs
  - *count_fastq_MB.sh* : counts the paired-end fastq files
  - *hmmscan-parser_MB.py* : parse the hmm files
------------

#### __Programs/scripts for metagenomic analysis__
- __Taxonomic profiling using metaphlan2__
  - *parse_metaphlan_MB.py* : The program will parse the metaphlan2 output and output file that can be directly used in R for downstream analysis.
- __Parser for CARD database__
  - *parse_card_MB.py*
  - *parse_aro_obo_MB.py*
------------

#### __Programs/scripts for isolates analysis__
- __Assembly of Isolate genomes__
  - *sequence_assembly_MB.py* : Performs assembly using Spades, IDBA-UD and Velvet and also runs the quast to evaluate genome assembly any finally generate summary for each samples
- __Characterization of Plasmid sequences from Spades Assembled genomes__ 
  - *plasmid_identify_MB.py* : 
- __Pangenome anaylsis__
