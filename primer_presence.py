import argparse
import logging
import sys
import os
import re
import subprocess
import difflib
from Bio import SeqIO

def main():

    # Argument Parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--seqs", help="Input text file of read sequences", required=True)
    parser.add_argument("-p", "--primers", help="Input text file of primers to look for", required=True)
    parser.add_argument("-v", "--verbose", help="Creates logging file with information for debugging", required=False, action='store_true')
    parser.add_argument("-r", "--reporting", help="Generates reports containing sequences", required=False, action='store_true')
    args = parser.parse_args()

    # Logging
    log_name = args.seqs + '.log'
    if args.verbose:
        logging.basicConfig(filename=log_name, encoding='utf-8', level=getattr(logging, "DEBUG", None))
    else:
        logging.basicConfig(filename=log_name, encoding='utf-8', level=getattr(logging, "INFO", None))
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(handler)
    logging.info('Logging started!')

    # Verify input files
    seq_file = args.seqs
    p_file = args.primers
    if not os.path.isfile(seq_file):
        logging.error('Error: read seqs file does not exist!')
        return
    if not os.path.isfile(p_file):
        logging.error('Error: primers file does not exist!')
        return

    # Parse the input read sequences
    logging.info('Parsing input reads..')
    total_reads = 0
    uniq_reads = 0
    reads = {}
    for read in open(seq_file).read().splitlines():
        total_reads += 1
        if reads.get(read):
            reads[read] += 1
        else:
            reads[read] = 1
            uniq_reads += 1
    logging.info('Total/unique = ' + str(total_reads) + "/" + str(uniq_reads))

    # Parse primer sequences
    logging.info('Parsing primers sequences..')
    primers = []
    for primer in open(p_file).read().splitlines():
        primers.append(primer)
    # Deduplicate primer seqs
    primers = list(dict.fromkeys(primers))
    logging.info('Total unique primer seqs = ' + str(len(primers)))
    

    # Search for primer hits within reads
    logging.info('Beginning primer search..')
    loc = []
    hit_dict = {}
    for primer in primers:
        hit_list = []
        for read in reads.keys():
            regex = re.search(primer,read)
            if regex:
                hit_list.append(read)
                loc.append(regex.span()[0])
        
        hit_dict[primer] = hit_list

        for found in hit_list:
            reads.pop(found)
    
    if args.reporting:
        logging.info('Generating reports..')
        p_read_file = open("primer_report.txt", "w")
        for key in hit_dict.keys():
            p_read_file.write(str(key) + "|" + str(hit_dict[key]) + "\n")
        p_read_file.close()
        logging.info("Primer Report Written..")
        np_read_file = open("nonprimed.seqs","w")
        for r in reads:
            np_read_file.write(r + "\n")
        np_read_file.close()
        logging.info("Non-Primed Sequences Written..")
    else:
        logging.info("Skipping report generation..")

    total_nonp = 0
    for seq in reads.keys():
        total_nonp += reads[seq]
    avg_pos = sum(loc) / len(loc)
    logging.info("Average primer start position: " + str(avg_pos))
    logging.info("Unique read sequences without primer match: " + str(len(reads)))
    logging.info("Total reads without primer match: " + str(total_nonp))

    # Cleanup logging
    if not args.verbose:
        os.remove("primer_presence.log")

if __name__ == "__main__":
    main()
    