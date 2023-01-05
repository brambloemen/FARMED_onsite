import sys, argparse, re, csv, time, logging, pysam, math
from datetime import datetime
from datetime import timedelta
from collections import defaultdict
import numpy as np
import plotly.express as px

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Get coverage and alignment statistics by time')
    parser.add_argument("bam", metavar="bam", help="bam file to process")
    parser.add_argument("-t", metavar="threads", dest="threads", type=int, help="Amount of threads to use")

    return parser


def preprocess(bamfile, threads):
    """
    Determine start and end time; duration of seq run
    Input: pysam.Alignmentfile object
    """

    fmt = '%Y-%m-%dT%H:%M:%SZ'
    start_time, end_time = None, None
    processedReads = defaultdict(list)

    bamfile = pysam.AlignmentFile(bamfile, "rb", threads=threads)
    template_lengths = {ref["SN"]:ref["LN"] for ref in bamfile.header.to_dict()["SQ"]}
    
    n_reads = 0
    for read in bamfile:

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

        # preprocess the reads, store useful variables in memory
        seq_time = seq_time - start_time
        seq_time = math.ceil(seq_time.total_seconds()/3600)
        
        processedReads[seq_time].append(read.to_string())
        n_reads += 1
    logging.info(f"number of reads: {n_reads}")
    duration = end_time - start_time
    duration = math.ceil(duration.total_seconds()/3600)

    return(start_time, duration, processedReads, template_lengths)

def process_reads(reads):
    # if more threads than reads, empty chunks will occur
    if len(reads)==0:
        return (None, None)
    imSummary = defaultdict(lambda: defaultdict(list))
    imCoverages = defaultdict(list)
    for read in reads:
        read = pysam.AlignedSegment.fromstring(read, samheader)
            
        # unmapped reads
        if read.reference_name == None:
            reference_name = "Unmapped"
            if not imSummary:
                imSummary = {reference_name: [0 ,0, 0, 0, 0, 0]}
            if reference_name not in imCoverages.keys():
                imSummary[reference_name] = [0 ,0, 0, 0, 0, 0]
                reference_length = 0
                    
                imSummary[reference_name][0] = reference_length
                imSummary[reference_name][1] += 1
                imSummary[reference_name][2] += read.query_length
                imSummary[reference_name][3] = 0
                imSummary[reference_name][4] = 0
                imSummary[reference_name][5] = 0
                continue

        reference_length = template_lengths[read.reference_name]

        # if first time a timepoint/reference is encounterd: populate dictionaries
        if not imSummary:
            # summary stats (cumulative over time): [ref_length, n_reads, total_readlength, matching_bases, coverage_breadth, coverage_depth]
            imSummary = {read.reference_name: [0 ,0, 0, 0, 0, 0]}
                    
        if read.reference_name not in imCoverages.keys():
            imCoverages[read.reference_name] =  [np.zeros(reference_length), np.zeros(reference_length)]

        if read.reference_name not in imSummary.keys():
            imSummary[read.reference_name] = [0 ,0, 0, 0, 0, 0]
                    
        # coverage
        imCoverages[read.reference_name][0][read.reference_start:read.reference_end-1] += 1

        # template identity
        templatePos = read.reference_start

        for i, c in enumerate(read.cigartuples):
            if i == 0 and c[0] == "S":
                continue
            if c[0] != "=":
                templatePos += c[1] # move template position to start of new cigar operation
                continue
            if c[0] == "M":
                match_end = templatePos + c[1] - 1
                imCoverages[read.reference_name][1][templatePos:match_end] += 1
                templatePos += c[1]

        imSummary[read.reference_name][0] = reference_length
        imSummary[read.reference_name][1] += 1
        imSummary[read.reference_name][2] += read.query_length
        imSummary[read.reference_name][3] += read.get_cigar_stats()[0][7]
    
    return(imSummary, imCoverages)
 
def parse_bam():
    """
    For each timepoint, calculate coverage, other stats.
    Next calculation cumulative stats over time
    """
    fmt = '%Y-%m-%dT%H:%M:%SZ'
    coverages = {}
    summary_stats = {}
    timepoints = [*range(1,duration + 1, 1)]

    for t in timepoints:
        summary_stats[t] = {}
        logging.info(f"Parsing reads sequenced in hour: {t}")
        
        reads = processedReads.pop(t)
            
        imSum, imCov = process_reads(reads)
        # in case there are empty chunks:
        if (imSum, imCov) == (None, None):
            continue
        for r in imCov.keys():
            # very first encounter
            if r not in coverages.keys():
                coverages[r] = imCov[r]
            else:
                # coverage breadth
                coverages[r][0] += imCov[r][0]
            if r not in summary_stats[t].keys():
                summary_stats[t][r] = [0] * 6
                summary_stats[t][r][0] = imSum[r][0]
                summary_stats[t][r][1] = imSum[r][1]
                summary_stats[t][r][2] = imSum[r][2]
                summary_stats[t][r][3] = imSum[r][3]

            # # coverage breadth
            # coverages[r] += imCov[r][0]

        for r in summary_stats[t].keys():
            summary_stats[t][r][4] = 100 * np.count_nonzero(coverages[r][0])/summary_stats[t][r][0]
            # coverage depth
            summary_stats[t][r][5] += summary_stats[t][r][2]/summary_stats[t][r][0]
        
    return(summary_stats)


if __name__ == '__main__':
    logging.info(f"Process arguments, load and preprocess reads")
    args = parse_arguments().parse_args()
    bamfile = args.bam
    threads = args.threads

    start_time, duration, processedReads, template_lengths = preprocess(bamfile, threads)
    samheader = pysam.AlignmentFile(bamfile, "rb", threads=threads).header
    logging.info(f"Start: {start_time}, duration (hours): {duration}")
    
    summary_stats = parse_bam()
    
    header = ["Time", "Template", "Template_length", "n_reads", "Total_readlength", "n_match_bases", "Coverage", "Depth"]
    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(header)

    for t, ref_stats in summary_stats.items():
        for ref, stats in ref_stats.items():
            writer.writerow([t, ref] + stats)

    logging.info(f"Done")

