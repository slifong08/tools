#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import os
import pandas as pd
import subprocess
import sys


###
#   arguments
###

arg_parser = argparse.ArgumentParser(description="intersect .bed with hg38 self chain")

arg_parser.add_argument(
    "-b", "--bedfile", help='bed file w/ full path')
arg_parser.add_argument(
    "-o", "--outpath", type=str, help='path to write outputs to')

args = arg_parser.parse_args()

# CONSTANTS
BED = args.bedfile
OUTPATH = args.outpath


# %% FUNCTIONS


def get_self_file():
    
    SELF_HG38 = "/dors/capra_lab/data/ucsc/hg38/self/hg38_repeats_self_coor.bed"
    
    return SELF_HG38


def make_files(path, sample_id):
    
    raw = os.path.join(path, f"{sample_id}_x_self-raw.bed")
    clean = os.path.join(path, f"{sample_id}_x_self-clean.bed")
    
    return raw, clean


def intersect_self(bedf, self_f, self_data_raw):
    """
    intersect bed file w/ self hg38 file. 
    
    bedtools intersection args: 
        -f .5, require 50% .bed element overlaps self element. 
        -wao
    """

    if os.path.exists(self_data_raw) is False:
        
        print("running self intersection")
        
        cmd = f"bedtools intersect -a {bedf} -b {self_f} -f 0.5 -wao > {self_data_raw}"
        subprocess.call(cmd, shell=True)

    else:
        print("ran self intersection already")


def format_df(self_data_raw, self_data_clean):
    
    """
    return cleaned self data. 
    
    input 
        self_data_raw (.bed) intersection of .bed file w/ self file. 

    output
        self_data_clean (.bed) dataframe with cleaned columns. 
        
    method
        1. load and name columns in dataframe (if the clean file has not already been written)
        2. add self column w/ default False
        3. format self overlap column, filling in non-self overlaps w/ zeros
        4. edit self column w/ True if overlapping self element is >10 bp
        5. keep specific columns
        6. write self_data_clean file. 
    """
    
    if os.path.exists(self_data_clean) is False or os.path.getsize(self_data_clean) == 0:
        
        # format the output dataframe
        #1 col names
        names = ["#chr", "start", "end", "bin", "#chrself", "startself", "endself", "selfoverlap"]

        # format dataframe
        df = pd.read_csv(self_data_raw, sep='\t',
                         header = None, 
                         names = names, 
                         low_memory=False)
        print(list(df))

        #2 add id for self

        df["self"] = False  # make col for counting self chain overlap.
        
        #3
        df["selfoverlap"] = df["selfoverlap"].fillna(0) # fill in non-overlapping regions
        
        #4
        df.loc[df["selfoverlap"].map(int) >= 10, "self"] = True  # Condition step
        
        #5
        group_cols = ["#chr", "start", "end", "bin", "self" ]

        newdf = df[group_cols].drop_duplicates()

        #6 write the dataframes to a file
        newdf.to_csv(self_data_clean, sep='\t', index=False)
    else:
        print("cleaned already")


#%%
def main(argv):

    SELF_HG38 = get_self_file()
    path, file = os.path.split(BED)
    sample_id = file.split(".bed")[0]
    
    
    SELF_DATA_RAW, SELF_DATA_CLEAN =  make_files(OUTPATH, sample_id)

    intersect_self(BED, SELF_HG38, SELF_DATA_RAW)
    format_df(SELF_DATA_RAW, SELF_DATA_CLEAN)

    
if __name__ == "__main__":
    main(sys.argv[1:])




