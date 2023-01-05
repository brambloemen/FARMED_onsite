import sys, argparse, re, csv, time, logging, pysam, math
from datetime import datetime
from datetime import timedelta
from collections import defaultdict
import numpy as np
import pandas as pd

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')


def parse_arguments():
    parser = argparse.ArgumentParser(description='Retrieve alignments from alignment file 1 in alignment file 2; \nPrint results to standard output')
    parser.add_argument("-t", metavar="taxa.bam", dest="bam_taxa", type=str, help="Filepath of bam with alignment to taxa database")
#   parser.add_argument("-r", dest="taxa_res", type=str, help="Filepath of res file of KMA alignmnent to taxa database")
    parser.add_argument("-a", metavar="AMR.bam", dest="AMR_bam", type=str, help="Filepath of bam with alignments to resistance database")
    parser.add_argument("--threads", type=int, help="Number of threads used to process SAM/BAM file", default=1)
    parser.add_argument("--tc_taxa", type=float, help="Minimum template coverage for a taxa template to be accepted", default=0.0)
    parser.add_argument("--tid_taxa", type=float, help="Minimum template identity for a taxa template to be accepted", default=0.0)
    parser.add_argument("--qid_taxa", type=float, help="Minimum query identity for a taxa template to be accepted", default=0.0)
    parser.add_argument("--tlen_taxa", type=float, help="Minimum template length for a taxa template to be accepted", default=0.0)
    parser.add_argument("--tdep_taxa", type=float, help="Minimum template depth for a taxa template to be accepted", default=1.0)
    parser.add_argument("--tid_arg", type=float, help="Minimum Percentage Identity to AMR gene template", default=97.0)
    parser.add_argument("--tdep_arg", type=float, help="Minimum depth to AMR gene template", default=1.0)
    parser.add_argument("--rlen", type=float, help="Minimum read length", default=0.0)

    return parser

def parse_taxa(taxa_bam, threads, tc_taxa, qid_taxa, tid_taxa, tlen_taxa, tdep_taxa):
    
    """
    Loop over alignment to taxa, filter out unwanted templates, determine start time and duration
    """
    
    # Filter for taxa -> insufficient taxa templates will not be considered for AMR linking
    fp_taxa_res = re.sub("\..am", ".res", taxa_bam)
    Filtered_taxa = pd.read_csv(fp_taxa_res, sep="\t")
    Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Coverage"] >= tc_taxa]
    Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Coverage"] >= qid_taxa]
    Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Identity"] >= tid_taxa]
    Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_length"] >= tlen_taxa]
    Filtered_taxa = Filtered_taxa[Filtered_taxa["Depth"] >= tdep_taxa]
    Filtered_taxa = Filtered_taxa["#Template"].tolist()


    fmt = '%Y-%m-%dT%H:%M:%SZ'
    start_time, end_time = None, None

    logging.info(f"Parsing alignment to taxa: {taxa_bam}")
    align_taxa = pysam.AlignmentFile(taxa_bam, "rb", threads=threads)
    align_taxa_refseq = {ref["SN"]:ref["LN"] for ref in align_taxa.header.to_dict()["SQ"]}
    align_taxa_n_lines = 0

    # align_taxa dictionary: {read:[Template1, readlength, matching bases, Template1 total reference length]}

    align_taxa_reads = defaultdict(lambda: [str, 0, 0, 0])
    for read in align_taxa:
        align_taxa_n_lines += 1
        # skip alignment if taxa template did not pass filter
        if read.reference_name not in Filtered_taxa:
            continue
        if read.query_length <= args.rlen:
            continue
        # store read identifiers
        n = read.query_name
        Template1 = read.reference_name
        l = read.query_length
        cigarEQ = read.get_cigar_stats()[0][7]
        Template1_l = align_taxa_refseq.get(Template1, None)
        if 100*cigarEQ/l < args.tid_taxa or 100*cigarEQ/l < args.qid_taxa:
            continue

        seq_time = re.search("start_time=\d+-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", read.query_name).group()
        seq_time = re.sub("start_time=", "", seq_time)
        seq_time = datetime.strptime(seq_time, fmt)

        if not start_time:
            start_time = seq_time
            end_time = seq_time
        elif seq_time < start_time:
            start_time = seq_time
        elif seq_time > end_time:
            end_time = seq_time

        seq_time = seq_time - start_time
        seq_time = math.ceil(seq_time.total_seconds()/3600)

        duration = end_time - start_time
        duration = math.ceil(duration.total_seconds()/3600)

        align_taxa_reads[n] = [seq_time, Template1, l, cigarEQ, Template1_l]

    logging.info(f"Done parsing alignment 1. Processed {align_taxa_n_lines} alignments")
    return (align_taxa_reads, start_time, duration)
    

def parse_AMR(amr_bam, processed_taxa_reads, threads, tid_arg, tdep_arg):
    """
    Loop over alignment to resistance database
    """

    # Filter for AMR genes -> insufficient AMR gene templates will not be considered for AMR linking
    fp_AMR_res = re.sub("\..am", ".res", amr_bam)
    Filtered_AMR = pd.read_csv(fp_AMR_res, sep="\t")
    Filtered_AMR = Filtered_AMR[Filtered_AMR["Template_Identity"] >= tid_arg]
    Filtered_AMR = Filtered_AMR[Filtered_AMR["Depth"] >= tdep_arg]
    Filtered_AMR = Filtered_AMR["#Template"].tolist()


    # summary dictionary. structure: {hour: {AMR gene: {Species: [#reads, total readlength, #exactly matched bases, readlength mapped to AMR, exactly matched bases to AMR, AMR template length]}}}
    alignment_links = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0, 0, 0])))
    logging.info(f"Parsing alignment to AMR: {amr_bam}")
    align_amr = pysam.AlignmentFile(amr_bam, "rb", threads=threads)
    align_amr_refseq = {ref["SN"]:ref["LN"] for ref in align_amr.header.to_dict()["SQ"]}
    align_amr_n_lines = 0

    for read in align_amr:
        align_amr_n_lines += 1
        # read identifiers
        n = read.query_name
        # alignment stats
        # not aligned to AMR --> skip
        if read.reference_name == None: 
            continue
        if read.query_length <= args.rlen:
            continue
        # skip reads matching to template not passing AMR filter
        if read.reference_name not in Filtered_AMR:
            continue
        else:
            Template2 = read.reference_name
        l = read.query_length
        cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)
        Template2_l = align_amr_refseq.get(Template2, None)
        if 100*(cigarEQ + read.get_cigar_stats()[0][8])/read.reference_length < args.tid_arg:
            continue

        if n in processed_taxa_reads:
            seq_time, Template1, Template1_readl, Template1_cigarEQ, Template1_l = processed_taxa_reads.pop(n) # pop matched reads from amr_reads: unmapped reads will remain
            alignment_links[seq_time][Template1][Template2][0] += 1
            alignment_links[seq_time][Template1][Template2][1] += Template1_readl # will be total length of all reads aligned to both template 1 and 2
            alignment_links[seq_time][Template1][Template2][2] += Template1_cigarEQ
            alignment_links[seq_time][Template1][Template2][3] += cigarEQ 
            alignment_links[seq_time][Template1][Template2][4] = Template1_l
            alignment_links[seq_time][Template1][Template2][5] = Template2_l
        else:
            alignment_links[seq_time]["Unmapped"][Template2][0] += 1
            alignment_links[seq_time]["Unmapped"][Template2][1] += l # will be total length of all reads aligned to template 2, with no match for temp
            alignment_links[seq_time]["Unmapped"][Template2][2] = None
            alignment_links[seq_time]["Unmapped"][Template2][3] += cigarEQ 
            alignment_links[seq_time]["Unmapped"][Template2][4] = None
            alignment_links[seq_time]["Unmapped"][Template2][5] = Template2_l

    logging.info(f"Done parsing alignment 2. Processed {align_amr_n_lines} alignments")

    return  alignment_links


if __name__ == '__main__':
    logging.info(f"Process arguments, load and preprocess reads")
    args = parse_arguments().parse_args()

    processed_taxa_reads, start_time, duration = parse_taxa(args.bam_taxa, args.threads, args.tc_taxa, args.qid_taxa, args.tid_taxa, args.tlen_taxa, args.tdep_taxa)
    logging.info(f"Start: {start_time}, duration (hours): {duration}")
    
    alignment_links = parse_AMR(args.AMR_bam, processed_taxa_reads, args.threads, args.tid_arg, args.tdep_arg)
    
    header = ["Time_hours", "Species", "ARG", "n_reads", "Total_readlength", "n_match_bases", "n_match_bases_ARG", "length_species_template", "length_ARG"]
    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(header)

    for t, arg_stats in alignment_links.items():
        for arg, spec_stats in arg_stats.items():
            for spec, stats in spec_stats.items():
                writer.writerow([t, spec, arg] + stats)

    logging.info(f"Done")

