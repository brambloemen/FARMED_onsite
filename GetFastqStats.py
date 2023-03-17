import sys
import gzip
import math

if len(sys.argv) < 2:
    print("Usage: python script.py input.fastq")
    sys.exit(1)

def mean_quality(qual_str):
    mq = -10 * math.log(sum([10**((ord(c) - 33)/-10) for c in qual_str]) / len(qual_str), 10)
    return mq

def parse_fastq(file):
    """
    Parse an opened fastq file
    """
    print("quals\tlengths")
    for i, line in enumerate(file):
        if i % 4 == 1:
            read = line.strip()
            read_length = len(read)
        elif i % 4 == 3:
            mean_qual = mean_quality(line.strip())
            # output to sys.stdout
            print(f"{mean_qual}\t{read_length}")

# Open the fastq file
if sys.argv[1].endswith(".gz"):
    with gzip.open(sys.argv[1], "rt") as fastq:
        parse_fastq(fastq)
else:
    with open(sys.argv[1], "r") as fastq:
        parse_fastq(fastq)
