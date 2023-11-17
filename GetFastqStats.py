import sys
import gzip
import math
import multiprocessing as mp
from tqdm import tqdm
import argparse
from itertools import islice

if len(sys.argv) < 2:
    print("Usage: python script.py input.fastq")
    sys.exit(1)

def mean_quality(qual_str):
    mq = -10 * math.log(sum([10**((ord(c) - 33)/-10) for c in qual_str]) / len(qual_str), 10)
    return mq

def process_chunk(chunk):
    return [
        (mean_quality(chunk[i].strip()), len(chunk[i].strip()))
        for i in range(3, len(chunk), 4)
    ]

#TODO: lazy evaluate so not entire file should be in memory
def parse_fastq_chunk(file, chunk_size=50000, num_threads=1):
    total_lines = sum(1 for _ in file)
    file.seek(0)  # Reset file pointer to the beginning

    chunks = [list(islice(file, chunk_size)) for _ in range(chunk_size)]
    results = []
    with mp.Pool(num_threads) as pool, tqdm(total=total_lines, unit='lines') as pbar:
        for result in pool.imap_unordered(process_chunk, chunks):
            results.extend(result)
            pbar.update(chunk_size)

    return results

def main():
    parser = argparse.ArgumentParser(description='Process a FASTQ file and calculate quality statistics.')
    parser.add_argument('input_file', metavar='input.fastq', type=str, help='Input FASTQ file')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for parallel processing')
    args = parser.parse_args()

    # Open the fastq file
    if args.input_file.endswith(".gz"):
        with gzip.open(sys.argv[1], "rt") as fastq:
            final_results = parse_fastq_chunk(fastq, num_threads=args.threads)
    else:
        with open(args.input_file, "r") as fastq:
            final_results = parse_fastq_chunk(fastq, num_threads=args.threads)
    
    # Print the final results
    print("quals\tlengths")
    for result in final_results:
        print(f"{result[0]}\t{result[1]}")

if __name__ == "__main__":
    main()