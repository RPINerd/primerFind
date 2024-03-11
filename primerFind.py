"""
    Primer Finder | RPINerd, 08/23/23

    Find and summarize the presence of given primers (or any subseq) within a fastq file
    Input is a plain text primers file with one sequence per line, as well as a fastq file
    to look up the sequences within
"""

import argparse
import logging
import os
import re
import sys
from datetime import datetime

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--seqs", help="Input text file of read sequences", required=True)
    parser.add_argument("-p", "--primers", help="Input text file of primers to look for", required=True)
    parser.add_argument(
        "-v",
        "--verbose",
        help="Creates logging file with information for debugging",
        required=False,
        action="store_true",
        default="store_false",
    )
    parser.add_argument(
        "-r", "--reporting", help="Generates reports containing sequences", required=False, action="store_true"
    )
    args = parser.parse_args()

    return args


def log_init(verbose=False):
    log_name = "primerFinder_{}.log".format(datetime.now().strftime("%y%m%d_%I%M%S"))
    if verbose:
        logging.basicConfig(filename=log_name, encoding="utf-8", level=getattr(logging, "DEBUG", None))
    else:
        logging.basicConfig(filename=log_name, encoding="utf-8", level=getattr(logging, "INFO", None))
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(handler)
    logging.info("Logging started!")


def fq_parse(file):
    logging.info("Parsing input reads..")

    records = {}
    total_reads = 0
    uniq_reads = 0
    for record in SeqIO.parse(file, "fastq"):
        total_reads += 1
        if records.get(str(record.seq)):
            records[str(record.seq)].append(record.id)
        else:
            records[str(record.seq)] = [record.id]
            uniq_reads += 1
    logging.info("Total/unique = " + str(total_reads) + "/" + str(uniq_reads))
    return records


def pr_parse(file):
    logging.info("Parsing primers sequences..")

    primers = []
    for primer in open(file).read().splitlines():
        primers.append(primer)

    # Deduplicate primer seqs
    logging.info("Submitted {} primer sequences..".format(len(primers)))
    primers = list(dict.fromkeys(primers))
    logging.info("Total of {} unique sequences.".format(len(primers)))

    return primers


def main(args):
    # Verify input files
    seq_file = args.seqs
    p_file = args.primers
    if not os.path.isfile(seq_file):
        logging.error("Error: read seqs file does not exist!")
        return
    if not os.path.isfile(p_file):
        logging.error("Error: primers file does not exist!")
        return

    # Parse the input read sequences
    reads = fq_parse(seq_file)

    # Parse primer sequences
    primers = pr_parse(p_file)

    # TODO prime candidate for threading/parallel work
    # Search for primer hits within reads
    logging.info("Beginning primer search..")
    loc = []
    hit_dict = {}
    for primer in primers:
        hit_list = []
        for read in reads.keys():
            regex = re.search(primer, read)
            if regex:
                hit_list.append(read)
                loc.append(regex.span()[0])

        hit_dict[primer] = hit_list

        for found in hit_list:
            reads.pop(found)

    if args.reporting:
        logging.info("Generating reports..")
        p_read_file = open("primer_report.txt", "w")
        for key in hit_dict.keys():
            p_read_file.write(str(key) + "|" + str(hit_dict[key]) + "\n")
        p_read_file.close()
        logging.info("Primer Report Written..")
        np_read_file = open("nonprimed.seqs", "w")
        for r in reads:
            np_read_file.write(r + "\n")
        np_read_file.close()
        logging.info("Non-Primed Sequences Written..")
    else:
        logging.info("Skipping report generation..")

    total_nonp = 0
    for seq in reads.keys():
        total_nonp += len(reads[seq])
    avg_pos = sum(loc) / len(loc)
    logging.info("Average primer start position: " + str(avg_pos))
    logging.info("Unique read sequences without primer match: " + str(len(reads)))
    logging.info("Total reads without primer match: " + str(total_nonp))

    # Cleanup logging
    if not args.verbose:
        os.remove("primer_presence.log")


if __name__ == "__main__":
    args = parse_args()
    log_init(args.verbose)
    main(args)
