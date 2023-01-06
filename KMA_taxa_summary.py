import pandas as pd
import numpy as np
import logging
import argparse
import re
import sys

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
parser = argparse.ArgumentParser(description='Summarize KMA result files for one or more experiments')
parser.add_argument("input", nargs='+', metavar="--input", help="Filepath(s) of KMA result file(s) (.res with corresponding .mapstat)")
args = parser.parse_args()

# dictionary to couple .res filepath to .mapstat filepath
filepaths = {fp : re.sub(".res", ".mapstat", fp) for fp in args.input}
# read .res and .mapstat files, merge, and add Experiment variable
def merge_kmares(resfile, mapstatfile):
    try:
        experiment=re.search("\w+.res", resfile).group()
        experiment_name = re.sub(".res", "", experiment)
        KMAres = pd.read_csv(resfile, sep='\t')
        KMAmapstat = pd.read_csv(mapstatfile, sep='\t', skiprows=6)
        KMA = pd.merge(KMAmapstat, KMAres, how='outer', left_on='# refSequence', right_on='#Template')
        KMA["Experiment"] = str(experiment_name) # add column with experiment from which the data originates
        return(KMA)
    except:
        raise Exception("Input should be a .res file, that has a corresponding .mapstat file, which is produced using the -e option of KMA.")

data = []
for res, mapstat in filepaths.items():
    print(res)
    pd.read_csv(res, sep='\t')
    KMA = merge_kmares(res, mapstat)
    data.append(KMA)

# concatenate all data
KMAresults = pd.concat(data)

# print summarized results to stdout in tsv format
KMAresults.to_csv(sys.stdout, sep="\t", index=False)