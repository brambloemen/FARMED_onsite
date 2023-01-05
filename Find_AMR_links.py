import pysam
import pandas as pd
import argparse
import re
import csv
from collections import defaultdict
import time
import logging
import sys
import plotly.express as px

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
parser = argparse.ArgumentParser(description='Retrieve alignments from alignment file 1 in alignment file 2; \nPrint results to standard output')
parser.add_argument("-t", metavar="taxa.bam", dest="bam_taxa", type=str, help="Filepath of bam with alignment to taxa database")
# parser.add_argument("-r", dest="taxa_res", type=str, help="Filepath of res file of KMA alignmnent to taxa database")
parser.add_argument("-a", metavar="AMR.bam", dest="AMR_bam", type=str, help="Filepath of bam with alignments to resistance database")
parser.add_argument("-o", "--output", type=str, default="ARGlinks")
parser.add_argument("--threads", type=int, help="Number of threads used to process SAM/BAM file", default=1)
parser.add_argument("--tc_taxa", type=float, help="Minimum template coverage for a taxa template to be accepted", default=0.0)
parser.add_argument("--tid_taxa", type=float, help="Minimum template identity for a taxa template to be accepted", default=0.0)
parser.add_argument("--qid_taxa", type=float, help="Minimum query identity for a taxa template to be accepted", default=0.0)
parser.add_argument("--tlen_taxa", type=float, help="Minimum template length for a taxa template to be accepted", default=0.0)
parser.add_argument("--tdep_taxa", type=float, help="Minimum template depth for a taxa template to be accepted", default=0.0)
parser.add_argument("--tid_arg", type=float, help="Minimum Percentage Identity to AMR gene template", default=0.0)
parser.add_argument("--tdep_arg", type=float, help="Minimum depth to AMR gene template", default=0.0)
parser.add_argument("--rlen", type=float, help="Minimum read length", default=0.0)
args = parser.parse_args()

# Filter for taxa -> insufficient taxa templates will not be considered for AMR linking
# this will reduce the amount of sequencing reads that will be parsed, speeding up the script
fp_taxa_res = re.sub("\..am", ".res", args.bam_taxa)
Filtered_taxa = pd.read_csv(fp_taxa_res, sep="\t")
Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Coverage"] >= args.tc_taxa]
Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Coverage"] >= args.qid_taxa]
Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_Identity"] >= args.tid_taxa]
Filtered_taxa = Filtered_taxa[Filtered_taxa["Template_length"] >= args.tlen_taxa]
Filtered_taxa = Filtered_taxa[Filtered_taxa["Depth"] >= args.tdep_taxa]
Filtered_taxa = Filtered_taxa["#Template"].tolist()

# Filter for AMR genes -> insufficient AMR gene templates will not be considered for AMR linking
fp_AMR_res = re.sub("\..am", ".res", args.AMR_bam)
Filtered_AMR = pd.read_csv(fp_AMR_res, sep="\t")
Filtered_AMR = Filtered_AMR[Filtered_AMR["Template_Identity"] >= args.tid_arg]
Filtered_AMR = Filtered_AMR[Filtered_AMR["Depth"] >= args.tdep_arg]
Filtered_AMR = Filtered_AMR["#Template"].tolist()


"""
Loop over alignment to taxa
"""
logging.info(f"Parsing alignment 1: {args.bam_taxa}")
align1 = pysam.AlignmentFile(args.bam_taxa, "rb", threads=args.threads)
align1_refseq = {ref["SN"]:ref["LN"] for ref in align1.header.to_dict()["SQ"]}
align1_n_lines = 0

# align1 dictionary: {read:[Template1, readlength, matching bases, Template1 total reference length]}

align1_reads = defaultdict(lambda: [str, 0, 0, 0])
for read in align1:
    align1_n_lines += 1
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
    Template1_l = align1_refseq.get(Template1, None)
    if 100*cigarEQ/l < args.tid_taxa or 100*cigarEQ/l < args.qid_taxa:
        continue

    align1_reads[n] = [Template1, l, cigarEQ, Template1_l]
    
logging.info(f"Done parsing alignment 1. Processed {align1_n_lines} alignments")


"""
Loop over alignment to resistance database
"""

# summary dictionary. structure: {AMR gene: {Species: [#reads, total readlength, #exactly matched bases, readlength mapped to AMR, exactly matched bases to AMR, AMR template length]}}
alignment_links = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0, 0, 0]))
logging.info(f"Parsing alignment 2: {args.AMR_bam}")
align2 = pysam.AlignmentFile(args.AMR_bam, "rb", threads=args.threads)
align2_refseq = {ref["SN"]:ref["LN"] for ref in align2.header.to_dict()["SQ"]}
align2_n_lines = 0

for read in align2:
    align2_n_lines += 1
    # read identifiers
    n = read.query_name
    # alignment stats
    # not aligned to AMR --> skip
    if read.reference_name == None: 
        continue
    # skip reads matching to template not passing AMR filter
    if read.reference_name not in Filtered_AMR:
        continue
    if read.query_length <= args.rlen:
        continue
    else:
        Template2 = read.reference_name
    l = read.query_length
    cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)
    Template2_l = align2_refseq.get(Template2, None)
    if 100*(cigarEQ + read.get_cigar_stats()[0][8])/read.reference_length < args.tid_arg:
        continue

    if n in align1_reads:
        Template1, Template1_readl, Template1_cigarEQ, Template1_l = align1_reads.pop(n) # pop matched reads from amr_reads: unmapped reads will remain
        alignment_links[Template1][Template2][0] += 1
        alignment_links[Template1][Template2][1] += Template1_readl # will be total length of all reads aligned to both template 1 and 2
        alignment_links[Template1][Template2][2] += Template1_cigarEQ
        alignment_links[Template1][Template2][3] += cigarEQ 
        alignment_links[Template1][Template2][4] = Template1_l
        alignment_links[Template1][Template2][5] = Template2_l
    else:
        alignment_links["Unmapped"][Template2][0] += 1
        alignment_links["Unmapped"][Template2][1] += l # will be total length of all reads aligned to template 2, with no match for temp
        alignment_links["Unmapped"][Template2][2] = None
        alignment_links["Unmapped"][Template2][3] += cigarEQ 
        alignment_links["Unmapped"][Template2][4] = None
        alignment_links["Unmapped"][Template2][5] = Template2_l

logging.info(f"Done parsing alignment 2. Processed {align2_n_lines} alignments")


# print results to standard output
logging.info(f"Writing results to tsv")
header = ["Template1", "Template2", "n_reads", "Total_readlength", "n_match_bases1", "n_match_bases2", "Template1_length", "Template2_length"]
writer = csv.writer(sys.stdout, delimiter="\t")
writer.writerow(header)
ARG_links = []

for temp1, temp2 in alignment_links.items():
    for template2, stats in temp2.items():
        if stats[0] < args.tdep_arg:
            continue
        ARG_links.append([temp1, template2] + stats)
        
        writer.writerow([temp1, template2] + stats)

logging.info(f"Done")

"""
Generate plots
"""
# extract species name
def extract_species_name(string):
    pattern = re.compile(r'[A-Z][a-z]+ [a-z]+')

    species = re.sub(r'\[|\]|\(|\)', '', string)
    species = re.sub(r'Synthetic|scf', '', species)
    if not re.search(pattern, species):
        species = "Unknown"
    else:
        species = re.search(pattern, species).group()

    return species

ARG_links = pd.DataFrame(ARG_links, columns=header)
ARG_links["Species"] = ARG_links["Template1"].apply(extract_species_name)
ARG_links["ARG"] = ARG_links["Template2"]
ARG_links = ARG_links.groupby(['ARG', 'Species'])['n_reads'].sum()
ARG_links = ARG_links.reset_index()
# Create the scatterplot
fig = px.scatter(ARG_links, x="Species", y="ARG", size="n_reads", size_max=50)

# Save the figure to a PNG file
fig.write_html(f"{args.output}_ARGlink.html")
fig.write_image(f"{args.output}_ARGlink.png")