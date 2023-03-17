import pandas as pd
import numpy as np
import logging
import argparse
import re
import sys
from Template_parsing import extract_species_name

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(message)s',
                    datefmt='%x %X')
def parser():
    parser = argparse.ArgumentParser(description='Summarize KMA result files for one or more experiments')
    parser.add_argument("input", nargs='+', metavar="--input", help="Filepath(s) of KMA result file(s) (.res with corresponding .mapstat)")
    args = parser.parse_args()
    # dictionary to couple .res filepath to .mapstat filepath
    filepaths = {fp : re.sub(".res", ".mapstat", fp) for fp in args.input}
    return(filepaths)
    
    
# read .res and .mapstat files, merge, add Experiment variable, and aggregate by species
def merge_kmares(resfile, mapstatfile):

    try:
        experiment=re.search("\w+.res", resfile).group()
        experiment_name = re.sub(".res", "", experiment)
        KMAres = pd.read_csv(resfile, sep='\t')
        KMAmapstat = pd.read_csv(mapstatfile, sep='\t', skiprows=6)
        KMA = pd.merge(KMAmapstat, KMAres, how='outer', left_on='# refSequence', right_on='#Template')
        KMA["Experiment"] = str(experiment_name) # add column with experiment from which the data originates
        KMA["Species"] = KMA["# refSequence"].apply(extract_species_name)
        KMA["Total_bp"] = sum(KMA["bpTotal"])
        KMA["Total_readCount"] = sum(KMA["readCount"])
        KMA = KMA.dropna(subset = ["bpTotal"])
        KMA = KMA[KMA['bpTotal'] > 0]
        grouped = KMA.groupby(['Species', 'Experiment'])

        # Calculate the weighted mean (by mapped bp) using the np.average function
        def nanaverage(A):
            dropped = A.dropna()
            return np.average(dropped, weights=KMA.loc[dropped.index, "bpTotal"])

        mean_queryID = grouped['Query_Identity'].agg(nanaverage)
        mean_query_coverage = grouped['Query_Coverage'].agg(nanaverage)
        mean_templateID = grouped['Template_Identity'].agg(nanaverage)
        mean_template_coverage = grouped['Template_Coverage'].agg(nanaverage)

        # Calculate other summary statistics
        p_bpTotal = grouped['bpTotal'].sum() / sum(KMA["bpTotal"])
        total_template_length = grouped['Template_length'].sum()
        mean_template_length = grouped['Template_length'].mean()
        depth = grouped['bpTotal'].sum() / mean_template_length
        p_readCount = grouped['readCount'].sum() / KMA['Total_readCount'].sum()
        readCount = grouped['readCount'].sum()
        mean_readlength = grouped['bpTotal'].sum() / grouped['readCount'].sum()
        refConsensusSum = grouped['refConsensusSum'].sum()

        # Combine the results into a single dataframe
        KMA = pd.DataFrame({'p_bpTotal': p_bpTotal,
                            'mean_queryID': mean_queryID,
                            'mean_query_coverage': mean_query_coverage,
                            'mean_templateID': mean_templateID,
                            'total_template_length': total_template_length,
                            'mean_template_length': mean_template_length,
                            'mean_template_coverage': mean_template_coverage,
                            'bpTotal': grouped['bpTotal'].sum(),
                            'depth': depth,
                            'p_readCount': p_readCount,
                            'readCount': readCount,
                            'mean_readlength': mean_readlength,
                            'refConsensusSum': refConsensusSum})

        KMA = KMA.reset_index()

        return(KMA)
    except:
        raise Exception("Input should be filepath to one or more .res files, with corresponding .mapstat files (produced using the -e option of KMA)")

if __name__ == '__main__':
    filepaths = parser()
    data = []
    # concatenate kma summary files from different experiments
    for res, mapstat in filepaths.items():
        # print(res)
        pd.read_csv(res, sep='\t')
        KMA = merge_kmares(res, mapstat)
        data.append(KMA)

    # concatenate all data
    KMAresults = pd.concat(data)

    # print summarized results to stdout in tsv format
    KMAresults.to_csv(sys.stdout, sep="\t", index=False)