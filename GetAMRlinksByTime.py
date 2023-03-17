import sys, argparse, re, csv, time, logging, pysam, math
from datetime import datetime
from datetime import timedelta
from collections import defaultdict
import numpy as np
import pandas as pd
import plotly.express as px
from Template_parsing import extract_species_name, extract_arg_name


logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')


def parse_arguments():
    parser = argparse.ArgumentParser(description='Retrieve alignments from alignment file 1 in alignment file 2; \nPrint results to standard output')
    parser.add_argument("-t", metavar="taxa.bam", dest="bam_taxa", type=str, help="Filepath of bam with alignment to taxa database")
#   parser.add_argument("-r", dest="taxa_res", type=str, help="Filepath of res file of KMA alignmnent to taxa database")
    parser.add_argument("-a", metavar="AMR.bam", dest="AMR_bam", type=str, help="Filepath of bam with alignments to resistance database")
    parser.add_argument("-o", metavar="output", dest="output", type=str, help="Name/filepath of output plots", default=None)
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


def extract_read_time(readname):
    fmt = '%Y-%m-%dT%H:%M:%SZ'
    seq_time = re.search("start_time=\d+-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", readname).group()
    seq_time = re.sub("start_time=", "", seq_time)
    seq_time = datetime.strptime(seq_time, fmt)
    return seq_time
    
def extract_exp_runtime(alignmentfile, threads):
    reads = pysam.AlignmentFile(alignmentfile, "rb", threads=threads)
    start_time, end_time = None, None
    for read in reads:
        seq_time = extract_read_time(read.query_name)
        if not start_time:
            start_time = seq_time
            end_time = seq_time
        elif seq_time < start_time:
            start_time = seq_time
        elif seq_time > end_time:
            end_time = seq_time
    return (start_time, end_time)


def parse_taxa(taxa_bam, threads, tc_taxa, qid_taxa, tid_taxa, tlen_taxa, tdep_taxa, start_time):
    
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

    logging.info(f"Parsing alignment to taxa: {taxa_bam}")
    align_taxa = pysam.AlignmentFile(taxa_bam, "rb", threads=threads)
    align_taxa_refseq = {ref["SN"]:ref["LN"] for ref in align_taxa.header.to_dict()["SQ"]}
    align_taxa_n_lines = 0

    # align_taxa dictionary: {read:[Template1, readlength, matching bases, Template1 total reference length]}
    align_taxa_reads = defaultdict(lambda: [str, 0, 0, 0])
    for read in align_taxa:
        align_taxa_n_lines += 1
        
        # skip if readlength too short
        if read.query_length <= args.rlen:
            continue
        # extract timing of read, calculate start and end time of sequencing experiment
        seq_time = extract_read_time(read.query_name)
        seq_time = seq_time - start_time
        seq_time = math.ceil(seq_time.total_seconds()/3600)

        # skip other steps if taxa template did not pass filter
        if read.reference_name not in Filtered_taxa:
            continue
        
        # store read identifiers
        n = read.query_name
        Template1 = read.reference_name
        l = read.query_length
        cigarEQ = read.get_cigar_stats()[0][7]
        Template1_l = align_taxa_refseq.get(Template1, None)
        if 100*cigarEQ/l < args.tid_taxa or 100*cigarEQ/l < args.qid_taxa:
            continue

        align_taxa_reads[n] = [seq_time, Template1, l, cigarEQ, Template1_l]

    logging.info(f"Done parsing alignment 1. Processed {align_taxa_n_lines} alignments")
    return align_taxa_reads
    

def parse_AMR(amr_bam, processed_taxa_reads, threads, tid_arg, tdep_arg, start_time):
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
        
        if read.query_length <= args.rlen:
            continue
        # read identifiers
        n = read.query_name
        # alignment stats
        # not aligned to AMR --> skip
        if read.reference_name == None: 
            continue
        
        # skip reads matching to template not passing AMR filter
        if read.reference_name not in Filtered_AMR:
            continue
        else:
            Template2 = read.reference_name
        l = read.query_length
        cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)
        Template2_l = align_amr_refseq.get(Template2, None)
        if 100*(cigarEQ)/read.reference_length < args.tid_arg:
            # + read.get_cigar_stats()[0][8]
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
            # read not in parsed reads => seq time not calculated
            seq_time = extract_read_time(n)
            seq_time = seq_time - start_time
            seq_time = math.ceil(seq_time.total_seconds()/3600)

            alignment_links[seq_time]["Unmapped"][Template2][0] += 1
            alignment_links[seq_time]["Unmapped"][Template2][1] += l # will be total length of all reads aligned to template 2, with no match for temp
            alignment_links[seq_time]["Unmapped"][Template2][2] = None
            alignment_links[seq_time]["Unmapped"][Template2][3] += cigarEQ 
            alignment_links[seq_time]["Unmapped"][Template2][4] = None
            alignment_links[seq_time]["Unmapped"][Template2][5] = Template2_l

    logging.info(f"Done parsing alignment 2. Processed {align_amr_n_lines} alignments")

    return  alignment_links

def create_output(output, alignment_links):

    header = ["Time_hours", "Species", "ARG", "n_reads", "Total_readlength", "n_match_bases1", "n_match_bases2", "Template1_length", "Template2_length"]
    # writer = csv.writer(sys.stdout, delimiter="\t")
    # writer.writerow(header)

    ARG_links = []

    for t, spec_stats in alignment_links.items():
        for spec, arg_stats in spec_stats.items():
            for arg, stats in arg_stats.items():
                # writer.writerow([t, spec, arg] + stats)
                ARG_links.append([t, spec, arg] + stats)

    """
    Generate plots
    """

    ARG_links = pd.DataFrame(ARG_links, columns=header)
    ARG_links["Species"] = ARG_links["Species"].apply(extract_species_name)
    ARG_links["ARG"] = ARG_links["ARG"].apply(extract_arg_name)
    
    grouped = ARG_links.groupby(['Time_hours', 'Species', 'ARG'])

    # Calculate other summary statistics
    n_reads = grouped['n_reads'].sum()
    Total_readlength = grouped['Total_readlength'].sum()
    n_match_bases1 = grouped['n_match_bases1'].sum()
    n_match_bases2 = grouped['n_match_bases2'].sum()
    Template1_length = grouped['Template1_length'].mean()
    Template2_length = grouped['Template2_length'].mean()
    Organism_QID = n_match_bases1/Total_readlength
    ARG_TID = n_match_bases2/(Template2_length * n_reads)

    ARG_links = pd.DataFrame({'n_reads': n_reads,
                        'Total_readlength': Total_readlength,
                        'n_match_bases1': n_match_bases1,
                        'n_match_bases2': n_match_bases2,
                        'Template1_length': Template1_length,
                        'Template2_length': Template2_length,
                        'Organism_QID': Organism_QID,
                        'ARG_TID': ARG_TID})

    ARG_links = ARG_links.reset_index()
    ARG_links.to_csv(sys.stdout, sep="\t", index=False)

    ARG_links["Species-ARG link"] =  ARG_links["Species"] + "-" + ARG_links["ARG"]
    ARG_links = ARG_links.groupby(["Time_hours", "Species-ARG link"])['n_reads'].sum()
    # ARG_links = ARG_links.sort_values(by='Time_hours')
    ARG_links = ARG_links.reset_index()
    ARG_links['n_reads'] = ARG_links.groupby(['Species-ARG link'])['n_reads'].cumsum()
    ARG_links = ARG_links.reset_index()

    # Create the plot
    fig = px.line(ARG_links, x='Time_hours', y='n_reads', color="Species-ARG link")

    # Save the figure to a PNG file
    if output != None:
        fig.write_html(f"{output}_ARGbytime.html")
        fig.write_image(f"{output}_ARGbytime.png")
    else:
        fig.write_html(f"ARGlink_by_time.html")
        fig.write_image(f"ARGlink_by_time.png")

if __name__ == '__main__':
    logging.info(f"Process arguments, load and preprocess reads")
    args = parse_arguments().parse_args()

    start_time, end_time = extract_exp_runtime(args.bam_taxa, args.threads)
    duration = end_time - start_time
    duration = math.ceil(duration.total_seconds()/3600)
    logging.info(f"Experiment start: {start_time}, duration (hours): {duration}")

    processed_taxa_reads = parse_taxa(args.bam_taxa, args.threads, args.tc_taxa, args.qid_taxa, args.tid_taxa, args.tlen_taxa, args.tdep_taxa, start_time)
        
    alignment_links = parse_AMR(args.AMR_bam, processed_taxa_reads, args.threads, args.tid_arg, args.tdep_arg, start_time)

    create_output(args.output, alignment_links)
    logging.info(f"Done")

