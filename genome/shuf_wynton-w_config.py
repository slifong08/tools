#!/usr/bin/env python
# coding: utf-8


from joblib import Parallel, delayed
import numpy as np
import os
import sys

sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")
sys.path.append("/wynton/home/ahituv/fongsl/tools/genome/")

import chr_functions
import config_readwrite as crw
import split_filename

CONFIG_NAME = sys.argv[1]  # full path to config file
SECTION = sys.argv[2]

config, configfile_name = crw.read_config(CONFIG_NAME)


"""
Output: N iterations of shuffled.bed with length-matched shuffled regions, chromosome-matched, non-overlapping, including regions specificed from incl arg

Input: .bed file, n shuffle iterations, genome build, background, sharedAcc, 50bp =< 0000 regions to include as the shuffle

Functions: 
- bedtools shuffle. 

Notes
- length match
- chromosome match
- include = only shuffle regions inside the included file. 


"""

###
#   arguments
###


TEST_BED = config[SECTION]["BED"]
ITERS =  config[SECTION]["ITERS"]
BUILD =  config[SECTION]["BUILD"]
INCLUDE =  config[SECTION]["INCL"]
SHUF_PATH =  config[SECTION]["OUTDIR"]

PATH, FILENAME, SAMPLE_ID = split_filename.split_filename(TEST_BED) 


###
#   functions
###


def loadConstants(build):  
    path_dict = {
                'hg19': ("/wynton/home/ahituv/fongsl/dna/black_list/hg19-blacklist.v2.bed", "/wynton/home/ahituv/fongsl/dna/hg19/hg19.chrom.sizes"),
                'hg38': ("/wynton/home/ahituv/fongsl/dna/black_list/hg38-blacklist.v2.gencode-exons.v42.bed", "/wynton/home/ahituv/fongsl/dna/hg38/hg38.chrom.sizes"), # exclude gencode exons, blacklist regions
                    }
    
    blacklist, sizes = path_dict[build]
    
    return blacklist, sizes


def shuffle(test_bed, shuf_path, sample_id, iter_, build, include):
    
    out = os.path.join(shuf_path, f"shuf-{sample_id}-{iter_}.bed")
    
    
    BLACKLIST, CHROM_SZ = loadConstants(build)  # not using BLACKLIST variable
        # instead, using include input as exclude if it has "blacklist"
        
    if "blacklist" in include:
        BEDshuf = f"bedtools shuffle -i {test_bed} -g {CHROM_SZ} -excl {include} -chrom -noOverlapping -maxTries 5000 > {out}" 
        
    else: 
        BEDshuf = f"bedtools shuffle -i {test_bed} -g {CHROM_SZ} -noOverlapping -maxTries 5000 -incl {include} > {out}"
        

    print(BEDshuf)
    
    os.system(BEDshuf)


###
#   Main
###


def main(argv):
    
    
    num_cores = 16
    print("number of cores", num_cores, "number of iters", ITERS)

    # run parallel jobs

    Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(shuffle)(TEST_BED, SHUF_PATH, SAMPLE_ID, i, BUILD, INCLUDE) for i in range(int(ITERS)))


    
if __name__ == "__main__":
    main(sys.argv[1:])