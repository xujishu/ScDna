import argparse
import pandas as pd
import numpy as np
import json
from os.path import basename
import sys


def merge_count_table(filelist,output):
    files = pd.read_csv(filelist,header=None)
    df = pd.DataFrame(files)
    merged = pd.DataFrame()
    for kk,ff in df.iterrows():
        filename = basename(ff[0])
        print(filename)
        sample_name = filename.split('_')[0]
        dat = pd.read_csv(
            ff[0],
            header=None,
            sep='\t',
            usecols=['ID','length',sample_name],
            names=['ID','length',sample_name,'unmapped'],
        )
        if kk == 0:
            merged = dat
        else:
            merged = pd.merge(left=merged,
                              right=dat[['ID',sample_name]],
                              left_on='ID',
                              right_on='ID'
                             )
    merged.to_csv(output, index=False)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--filelist",
        dest="filelist",
        required=True,
        help=
        "a file which list all input files to parse"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help=
        "a file to output"
    ) 
    args = parser.parse_args()
    merge_count_table(args.filelist,args.output)

if __name__ == "__main__":
            main()
