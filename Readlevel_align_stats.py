import pysam
import argparse
import re
import csv
from collections import defaultdict
import logging
from Template_parsing import extract_species_name

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
parser = argparse.ArgumentParser(description='Extract read-level and template-level mapping statistics from KMA-analyes, including unmapped reads\n (which are often not returned by KMA).')
parser.add_argument("bam", metavar="--input.bam", help="Filepath of bamfile", nargs="+")

args = parser.parse_args()

for file in args.bam:
    readstatsfp = re.sub(".bam|.sam", "", file) + ".readstats.class.tsv"
    with open(readstatsfp, "w") as csv_reads:

        csv_readstats = csv.writer(csv_reads, delimiter="\t")
        csv_readstats.writerow(["Read_ID", "Template", "Length", "Matched_bp"])

        # structure: [number of reads, combined length of reads, total sequence match bp]
        species_sum = defaultdict(lambda: [0, 0, 0])

        """
        Loop over samfile. Write required alignment info to csv files.
        """
        logging.info("Parsing BAM file")
        samfile = pysam.AlignmentFile(file, "rb", threads=1)
        samlines = 0
        for read in samfile:
            # Most simple Read ID
            n = str(read.query_name).split(" ")[0]
            n = re.sub("\D", "", n)

            # alignment stats
            s = extract_species_name(read.reference_name)
            l = read.query_length
            # cigarM += read.get_cigar_stats()[0][0] # alignement match (a base in query aligns with a base in template)
            cigarEQ = read.get_cigar_stats()[0][7] # EQ: exact sequence match (base in query matches base in template)

            # summarize stats by species
            species_sum[s][0] += 1
            species_sum[s][1] += l
            species_sum[s][2] += cigarEQ

            csv_readstats.writerow([n, s, l, cigarEQ])
            samlines += 1
