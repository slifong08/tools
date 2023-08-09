#!/usr/bin/env python
# coding: utf-8


import argparse
from collections import Counter
import configparser
import glob
import numpy as np
import os
import pandas as pd
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")

import config_readwrite as crw
import chr_functions


"""
Output: .tsv with sequence id,  GC count, and GC density

Input: .bed file.

Functions: 
- phast suite split_msa function

Notes


"""

###
#   arguments
###

arg_parser = argparse.ArgumentParser(description= "compute sequence identity between two species")

arg_parser.add_argument("-b","--bedfile", help='.bed file in species 1 coordinates')
arg_parser.add_argument("-s","--subject", help='subject - name species 1 (e.g. hg38)')
arg_parser.add_argument("-q","--query", help='query - name species 2 (e.g. rheMac10)')
arg_parser.add_argument("-o","--outdirectory", help='out directory to save results')
arg_parser.add_argument("-a","--alignment", help='options = hg38.rheMac10, multiz100way, multiz30way, multiz20way')

args = arg_parser.parse_args()

TEST_BED = args.bedfile
SUBJECT= args.subject
QUERY = args.query
OUTDIR = args.outdirectory 
ALIGNMENT = args.alignment


# write these directories, files

sample_id = os.path.splitext(os.path.basename(TEST_BED))[0]
chr_path = os.path.join(OUTDIR,  f"chr-{sample_id}")

SEQID_DATA_RAW = os.path.join(OUTDIR, f"{sample_id}-seq_identity_raw.tsv")
SEQONLY_DATA = os.path.join(OUTDIR, f"{sample_id}-seq_only.tsv")
SEQID_DATA = os.path.join(OUTDIR, f"{sample_id}-seq_identity.tsv")

PATHS = [OUTDIR, chr_path]
DATA_PATHS = ",".join(PATHS)

print(chr_path, sample_id, DATA_PATHS)

### 
# Functions
### 

def make_paths(data_paths):

    if type(data_paths) is list:
        data_paths = ",".join(data_paths)
        
    for i, path in enumerate(data_paths.split(",")):
 
        if os.path.exists(path) is False:

            os.mkdir(path)
            print("made", i, path)

def split_chr(test_bed, chr_path):

    chr_functions.split_into_chr_bed(test_bed, chr_path)

    
def get_maf_src(chr_, subject):
    
    # return species' chromosome.maf file
    
    maf_path = f'/dors/capra_lab/data/ucsc/{subject}'
    maf_chr = os.path.join(maf_path, ALIGNMENT, f"{chr_}.maf") 
        
    return maf_chr


def msa_split_x_bed(chr_, chr_bed, path, chr_raw_path): # msa
    
    """
    1. go to path of bed file
    2. get chr-specific maf file in subject genome
    3. write feat-arg string
    4. write outroot arg string
    5. write msa_split cmd string
    6. check that .fa files have not been split already. 
    7. if not, split
    8. return list of .fa files - there are N .fa files for N rows of the bed file. 
    split msa by bed file. 
    use phast Suite's msa_split function w/ --for-features to split on bedfile .
    
    """
    
    msa_split_bin = "/dors/capra_lab/bin/./msa_split"

    #1
    os.chdir(path)
    
    #2
    maf_file = get_maf_src(chr_, SUBJECT)
    
    #3
    feat_arg = f"--features {chr_bed} --for-features"  # --for-features will split on each row of bed file. 
    
    #4
    out_root_arg = f"--out-root {chr_}"
    
    #5
    cmd = f"{msa_split_bin} {maf_file} --in-format MAF {feat_arg} {out_root_arg}"
    
    #6
    already_split = len(glob.glob(f"{chr_}*.fa"))  # outputs lines (as individual files) 
   
    n_lines = sum(1 for line in open(chr_bed))  # input lines (from one file) 
    
    #7
    if already_split !=n_lines:
        print(cmd)
        subprocess.call(cmd, shell = True)
        
    else:
        print("done")
    #8
    already_split = glob.glob(f"{chr_}*.fa")
    
    return already_split



# msasplit indexes by 1. Need to reindex at zero. 
def reset_zero_index(coordinate):
    coordinate = int(coordinate) -1
    return coordinate



def get_percent_identity(subjSeq, querySeq):

    lenSeq = len(subjSeq) # get the length of the sequence alignment.

    count_identical = 0
    count_gap = 0
    count_non_identical = 0
    
    # parse through sequence and ask if alignments match. 
    for a,b in zip(subjSeq,querySeq):

        if a==b:
            count_identical+=1  # count identical bases

        elif a != b:
            count_non_identical +=1  # count non-identical bases
            
        if a == "-" or b == "-":
            count_gap +=1  # count gap bases
            
    percent = count_identical/lenSeq  # return percent identity

    return count_identical, count_gap, percent


# In[10]:


def make_region_df(fa_handle, sub_seq, qry_seq, subsize, qrysize, score, gap, percent):

    chr_, coor = fa_handle.split(".")[0:2]
    start, end = coor.split("-")

    start= reset_zero_index(start) #, reset_zero_index(end)  # 0-index instead of 1-index

    df = pd.DataFrame({  # make a dataframe of the results
    "#chr" :[chr_],
    "start":[start],
    "end":[end],
    f"{SUBJECT}_Seqlen": [subsize],
    f"{QUERY}_Seqlen": [qrysize],
    f"{SUBJECT}_seq":sub_seq,
    f"{QUERY}_seq":qry_seq,
    "score":[score],
    "gap":[gap],
    "percent_identity":[percent],

    })

    return df


# In[11]:


def make_block_df(results_dict, chr_, path):

    # concat the dictionary
    re = pd.concat(results_dict.values()).drop_duplicates()
    re_seq = re[['#chr', 'start', 'end',f'{SUBJECT}_seq', f'{QUERY}_seq']].copy()
    re = re[[
            '#chr', 'start', 'end',
            f'{SUBJECT}_Seqlen', f'{QUERY}_Seqlen',
            'score', 'gap', 'percent_identity',
            ]]
    
    # write identity outfile
    outf = f"{chr_}_seq_identity.tsv"
    out = os.path.join(path, outf)
    
    # write identity outfile
    outf_seq = f"{chr_}_seq_only.tsv"
    outseq = os.path.join(path, outf_seq)

    # write the files
    re.to_csv(out, sep = '\t', index = False)
    re_seq.to_csv(outseq, sep = '\t', index = False)

    return re


# In[12]:


def concat(path):
    
    """
    (1) go to path
    (2) concat chr*_seq_identity.tsv files
    (3) remove chromosome-specific seq_identity.tsv files
    (4) concat chr*_seq_only.tsv files
    (5) remove chromosome-specific seq_only.tsv files
    (6) intersect seq_id w/ original bed file
    """
    #(1)
    os.chdir(path)
    
    #(2)
    cmd = f'cat *_seq_identity.tsv > {SEQID_DATA_RAW}'
    
    #(3)
    if os.path.exists(SEQID_DATA_RAW) is False or os.path.getsize(SEQID_DATA_RAW) ==0:
        subprocess.call(cmd, shell = True)
    
        cmd = 'rm *_seq_identity.tsv' # clean up
        subprocess.call(cmd, shell = True)
    else:
        print("\nmade seqid_data\n")
    #(4) 
    cmd = f'cat *_seq_only.tsv > {SEQONLY_DATA}'
    #(5)
    if os.path.exists(SEQONLY_DATA) is False or os.path.getsize(SEQONLY_DATA) ==0:
        subprocess.call(cmd, shell = True)
    
        cmd = 'rm *_seq_only.tsv' # clean up
        subprocess.call(cmd, shell = True)
    else:
        print("\nmade seqonly_data\n")
    
    # (6) keep only the files that overlap 90% of an identity tile
    cmd = f"bedtools intersect -a {TEST_BED} -b {SEQID_DATA_RAW} -f 0.9 -wao > {SEQID_DATA}"
    subprocess.call(cmd, shell = True)
    

def extract_fa_data(fa_handle):
    if os.path.exists(fa_handle) is True:  # check that the path exists
        with open(fa_handle, "r") as fa_reader:
            """
            (1) set empty values for collecting species' sequence and sequence size
            (2) if species is human, set species variable to human
            (3) else, set species variable to rhesus
            (4) if neither hg38 or rheMac10 annotation, use species variable to recore sequence, size
            """
            #(1)

            sub_seq, qry_seq = "", ""
            subsize, qrysize = 0,0
            species = None

            for i, line in enumerate(fa_reader):


                #(2)
                if SUBJECT in line:
                    species = SUBJECT
                #(3)
                elif QUERY in line:
                    species = QUERY
                #(4)
                else:

                    line = line.strip("\n")  # strip the \n
                    if species == SUBJECT:
                        sub_seq += line
                        subsize += len(line)
                    elif species == QUERY:
                        qry_seq += line
                        qrysize += len(line)
                #os.remove(fa_handle)  # delete the handle

    else:
        print("no fa", fa_handle)
        sub_seq, qry_seq, subsize, qrysize = None, None, -1, -1
    return sub_seq, qry_seq, subsize, qrysize


def make_chr_files(chr_, working_dir, chr_path):
    outf = os.path.join(working_dir, f"{chr_}_seq_identity.tsv")  # chr-seq identity file to write
    chrF = f"{chr_}.bed"  # regions to do msasplit on. 
    chr_bed = os.path.join(chr_path, chrF)  # with full path. 

    return outf, chrF, chr_bed
 

def main(argv):

    chrList = chr_functions.make_chr_list()  # get chromosomes

    if os.path.exists(SEQID_DATA) is False:   

        """
        (0) Make output paths
        (1) split file by chromosome number
        """

        #(0)
        make_paths(DATA_PATHS)

        #(1)
        split_chr(TEST_BED, chr_path) 


        for chr_ in chrList:  
            """
            per chromosome this section will
            (2) perform msa_splits 
            (3) quantify sequence identity
            """

            print(chr_)
            
            outf, chrF, CHR_BED = make_chr_files(chr_, OUTDIR, chr_path)
            

            #(2) perform msa_splits
            if os.path.exists(outf) is False or os.path.getsize(outf) ==0:

                msa_splits = msa_split_x_bed(chr_, CHR_BED, OUTDIR, chr_path) # MSA
                print(msa_splits, "\n\n finished splitting \n\n")

                results_dict = {}  # for collecting data on sequence identity

                #(3) quantify sequence identity
                for n, fa_handle in enumerate(msa_splits):
                
                    if n%10000 == 0 and n !=0:
                        print(n)
                    
                    os.chdir(OUTDIR)
                    
                    sub_seq, qry_seq, subsize, qrysize = extract_fa_data(fa_handle)

                    os.remove(fa_handle) ## remove the handle
                    
                    count_identical, count_gap, percent = get_percent_identity(sub_seq, qry_seq)

                    df = make_region_df(fa_handle, sub_seq, qry_seq, subsize, qrysize, count_identical, count_gap, percent)
                    
                    results_dict[n] = df

                re = make_block_df(results_dict, chr_, OUTDIR)  # write the results ot a file


        """
        (4) consolidate each chromosomes' dataframe
        """
        re = make_block_df(results_dict, chr_, OUTDIR)

        """
        (5) concatenate all chromosome results and intersect w/ original file. 

        """
        concat(OUTDIR) 

if __name__ == "__main__":
    main(sys.argv[1:])