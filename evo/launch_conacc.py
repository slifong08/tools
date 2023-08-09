#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import argparse
import configparser
import glob
import itertools 
import numpy as np
import os, sys
import subprocess
import time
sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import chr_functions

###
# args
###

arg_parser = argparse.ArgumentParser(description="estimate conservation/acceleration w phyloP")

arg_parser.add_argument("-b", "--bedfile", help="bed file")

arg_parser.add_argument("-br", "--branches", help='hg38, rheMac3')
arg_parser.add_argument("-c", "--chromosome_directory", help='hg38, rheMac3')

arg_parser.add_argument(
    "-msa", "--multiz", help='20-, 30-, 100-way multiz in hg38')

arg_parser.add_argument(
    "-mod", "--model", help='full", hg38-rheMac8', default="full")
arg_parser.add_argument("-o", "--outdirectory", help = "outdirectory to dump results")



# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

BEDF = args.bedfile
MSAWAY = args.msa
GENOME_BUILD = args.msa
DATA_PATH = args.outdirectory
CHR_PATH = args.chromosome_directory
# CONSTANTS
ARRAY = True
BIN_PATH = '/dors/capra_lab/users/fongsl/tools/evo/'    

##
# FUNCTIONS
###

def get_phylop_dict(msaway, genome_build):

    phylop_dict = {
        'phylop_bin':'/dors/capra_lab/bin/./phyloP', 
        'dors_maf_path': f'/dors/capra_lab/data/ucsc/{genome_build}',
        'maf':f'/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/maf/{msaway}.maf.gz', 
        'branches':['hg38', 'rheMac8', 'hg38-rheMac8'], 
        'models': ['full', 'rheMac8_noOWM', 'hg38_noAPES'],
    }
    
    neutral_model_dict = {
                            'full': f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}.mod', 
                            'hg38-rhemac8': f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}_hg38-rheMac8.mod',
                            'rhemac8_noowm': f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}_rheMac8_noOWM.mod', 
                            'hg38_noapes':f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}_hg38_noAPES.mod'
                            }
    return phylop_dict, neutral_model_dict


# in case you need to split file on size before getting started
def split_by_line(f, data_path, chr_num):

    chr_path = os.path.join(data_path, chr_num) # make dir for chromosome splits
    
    try:
        os.mkdir(chr_path)
    except FileExistsError:
        pass
    
    # change dir to the output chr path (not the original CHR_PATH, 
    # where file is split on CHR, but not line number)
    os.chdir(chr_path)

    small_fs = glob.glob(f"{chr_path}/{chr_num}-*")

    # split the file in command line into sizes of 1000 lines
    cmd = f"split -l 1000 {f} {chr_num}-"
    if len(small_fs) ==0:
        print("splitting")
        subprocess.call(cmd, shell = True)

    else:
        print("already split")
    small_fs = glob.glob(f"{chr_path}/{chr_num}-*")

    return small_fs


def make_run_list(branches, models, chrs):
    runs = []
    no_runs = [('hg38', 'rheMac8_noOWM') , 
               ('hg38-rheMac8', 'rheMac8_noOWM'), 
               ('hg38-rheMac8', 'hg38_noAPES'),
              ('rheMac8', 'hg38_noAPES') 
              ] # don't run these tuples. Not interested yet in these results
    for b in branches:
        for m in models:
            for c in chrs:           
                combo = [b, m, c]
                if combo not in runs and (b,m) not in no_runs:
                    runs.append(combo)
                    
    return runs

def print_cmd(bin_path, branches, msaway, mod, chrnum, jobtype, file):

    # tell us what is being run
    print("\nrunning", jobtype,  "\non", file,
    "\nbranches:", branches, "\nmsa:", msaway,
    "\nmod:", mod, "\nchr:", chrnum, "\n")

    
def run_conacc_slurm_array(bin_path, chrnum, branches, msaway, mod, data_path, genome_build):

    script = os.path.join(bin_path, "conacc_array.slurm")
    
    chr_path = os.path.join(data_path, chrnum)
    
    num_files = len(os.listdir(os.path.join(data_path, chrnum)))# get the number of files and set the array arg. 
    

    big_mem_chr = ["chr1", "chr2", "chr6", "chr6", "chr5", "chr11", "chr7", "chr12","chr17","chr19",]

    if chrnum in big_mem_chr:
        mem = "--mem=120GB"
        array = f"--array [0-{num_files}:5]%5"
        print("mem requested", mem)
    else:
        mem = "--mem=64GB"
        array = f"--array [0-{num_files}:10]%10"

    # make the command

    cmd = f"sbatch {array} {mem} {script} {chrnum} {genome_build} {branches} {msaway} {mod} {data_path}"
    print(cmd)
    jobtype = "slurm"
    file = "array"
    print_cmd(bin_path, branches, msaway, mod, chrnum, jobtype, file)
    # run it
    return cmd


# function to check if you have already run these files. 

def check_already_run(combo, msaway, phylop_path, chr_path):
    
    branch, model, chrnum = combo[0], combo[1], combo[2].split(".bed")[0]
      
    outpath = os.path.join(
    phylop_path, f"{chrnum}", f"multiz{msaway}way_br-{branch}_mod-{model}")

    if os.path.exists(outpath) is True:  # if you have tried to run this before. 

        
        finished = os.path.join(outpath, f"{chrnum}_conacc.bed")  # if the run was complete, there should be a file named this. 
        
        if os.path.exists(finished) is True:  # now, make sure that the file is exactly the size it should be. 
           
            
            finished_lines = sum(1 for line in open(finished))
            expected_lines = sum(1 for line in open(os.path.join(chr_path, f"{chrnum}.bed")))
            
            if finished_lines == expected_lines:  # if obs result lines is the expected number of lines.  
                
                skip_run = True
            else:
                skip_run = False
                print("finished v. expected", finished_lines, expected_lines)
        else:
            print("never finished", finished)
            skip_run = False
    else:
        print("never made folder and never run", outpath)
        skip_run = False
    return skip_run



def main(argv):
    
    PHYLOP_DICT, NEUTRAl_MODELS = get_phylop_dict(MSAWAY, GENOME_BUILD)
    BRANCHES = PHYLOP_DICT["branches"] # get from dict
    MODELS = PHYLOP_DICT["models"]

    os.chdir(CHR_PATH)
    chrs = glob.glob("chr*.bed")  # get all the chromosome files

    # exclude these chromosomes by taking set difference.
    excl_chr = set(['chrX.bed', 'chrY.bed', 'chrM.bed', 
                    'chr14_KI270726v1_random.bed', 'chr16_KI270728v1_random.bed'
                    'chr14_KI270722v1_random.bed', 'chr16_KI270728v1_random.bed', 
                    'chr14_KI270722v1_random.bed'
                   ])
    
    chrs_ = list(set(chrs).difference(excl_chr))

    
    for chr_ in chrs_: 

        chr_ = chr_.split(".bed")[0]

        chr_f = os.path.join(CHR_PATH, (chr_+".bed"))
        
        # split chr_files into 1000 line files.
        split_fs = split_by_line(chr_f, DATA_PATH, chr_)

    run_list = make_run_list(BRANCHES, MODELS, chrs_)
    
    print(len(run_list),"phylop chr x branch x model comparisons to run", )

    val = 0

    # per chr-branch-model combination in run list

    for run in run_list[val:]:
        
        BRANCH = run[0]
        MODEL = run[1]
        CHR = run[2]
       
        SKIP_RUN = check_already_run(run, MSAWAY, DATA_PATH, CHR_PATH)  # check if you have already run this:

            
        cmd = run_conacc_slurm_array(BIN_PATH, os.path.join(CHR_PATH, CHR), GENOME_BUILD, BRANCH, MSAWAY, MODEL, DATA_PATH)

        """
        if (val % 9) == 0 and val>0:

            sleeptime = 60*25
            time.sleep(sleeptime)
        """    
                
        if SKIP_RUN is False:
            print(val, cmd)
            subprocess.call(cmd, shell = True)
        
            val +=1
        else:
            print("run already")

if __name__ == "__main__":
    main(sys.argv[1:])




