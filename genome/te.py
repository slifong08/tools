#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import configparser
import os
import pandas as pd
import subprocess
import sys
sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import split_filename


"""
Output: .tsv with repeatmasker te name, te fam results, summarized by bed region

Input: .bed file, genome build, outdirectory path

Functions: 
- bedtools intersect w/ repeatmasker coordinates + dataframe cleanup

Notes

- includes only TEs that overlap >= 10bp


"""

###
#   arguments
###

arg_parser = argparse.ArgumentParser(description= "compute sequence identity between two species")

arg_parser.add_argument("-b","--bedfile", help='.bed file in species 1 coordinates')
arg_parser.add_argument("-g","--genome_build", help='repeatmasker genome build (e.g. hg38)')
arg_parser.add_argument("-o","--outdirectory", help='out directory to save results')


args = arg_parser.parse_args()

BEDF = args.bedfile
BUILD = args.genome_build
OUTDIR = args.outdirectory


###
#   FUNCTIONS
###

def repeatmasker_files(build):
    repeatmasker_dict = {"dors_path":"/dors/capra_lab/data/transposable_elements/repeatmasker/",
                        "hg38":"/dors/capra_lab/data/transposable_elements/repeatmasker/hg38.bed",
                        "rheMac10":"/dors/capra_lab/data/transposable_elements/repeatmasker/rheMac10.bed"}

    return repeatmasker_dict[build]
    
# In[4]:


def intersect_te(bedf, te, build, data_path):
    path, filename, sample_id = split_filename.split_filename(bedf)

    print("Intersecting .bed x repeatmasker for", build)
    os.chdir(path)

    # get info to write new files

    out_te = os.path.join(data_path, f"{sample_id}_TE-raw.bed")

    # if you haven't done this intersection already...
    
    cmd = f"bedtools intersect -a {bedf} -b {te} -wao > {out_te}"    
    
    if os.path.exists(out_te) is False:
        print(cmd)
        subprocess.call(cmd, shell=True)
    else:
        print("already intersected w/ TE")

    return out_te, sample_id


def format_outdf(out_te, build, outdir, sample_id):
    """
    Summary: Format the output file, count TE overlaps, record TE names 
    
    1. name columns
    2. load dataframe
    3. count TE overlaps w/ >=10bp
          Condition: if TE overlaps window by 10bp or more,
          set 'te_count' column equal to one, 
          Count: groupby and sum TE counts
        
    4. groupby and concat on TE names, TE fam
    5. save file
    """
    print("formatting the TE intersection file ")

    #1
    if "rheMac10" in out_te:
        # col names
        cols = ["#chr", "start", "end", "bin",
                "chr_te", "start_te", "end_te", "te", 
                "te_fam", "strand", "len_te_overlap"]
        
    elif "Fong_regions.rheMac10_unmapped_TE" in out_te:
        cols = ["#chr", "start", "end","bin", 
                "chr_te", "start_te", "end_te", "te", 
                "strand", "len_te_overlap"]
   
    else:
        # col names
        cols = ["#chr", "start", "end", 'region_id',
                "chr_te", "start_te", "end_te", 
                "te", "te_fam", "len_te_overlap"]

    #2 load dataframe
    df = pd.read_csv(out_te, sep='\t', header=None, low_memory=False)
    print(out_te, df.head())

    # rename columns
    df.columns = cols
   
    #3 add id for te count
    
    idname = f"te_count-{build}"

    df[idname] = 0  # make col for counting overlapping TEs.

    df.loc[df["len_te_overlap"].map(int) >= 10, idname] = 1  # Condition step

    group_cols = list(df.columns[:4])

    newdf = df.groupby(group_cols)[idname].sum().reset_index()  # Count step

    #4 TE names step
    te_namedf = df.groupby(group_cols)["te"].unique().reset_index()
    
    #print("newdf", list(newdf))
    #print("te_namedf", list(te_namedf), te_namedf)
    
    te_fam = df.groupby(group_cols)["te_fam"].unique().reset_index()

    # merge the dataframes
    newdf = pd.merge(newdf, te_namedf, how="left")
    newdf = pd.merge(newdf, te_fam, how="left")

    #5 write the dataframes to a file
    outf = os.path.join(outdir, f"{sample_id}_TE-cleaned.bed")
    print(outf)

    newdf.to_csv(outf, sep='\t', index=False)

    return newdf


# In[5]:


def main(argv):  # %% Run the functions
    
    TE = repeatmasker_files(BUILD) # get the repeatmasker file to intersect
    OUTTE, SAMPLE_ID = intersect_te(BEDF, TE, BUILD, OUTDIR) # do the bedtools intersection
    
    OUTDF = format_outdf(OUTTE, BUILD, OUTDIR, SAMPLE_ID)  # format the dataframe

if __name__ == "__main__":
    main(sys.argv[1:])





