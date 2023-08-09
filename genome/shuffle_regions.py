#!/usr/bin/env python
# coding: utf-8


import argparse
from joblib import Parallel, delayed
import numpy as np
import os
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")

import config_readwrite as crw
import chr_functions
import split_filename


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

arg_parser = argparse.ArgumentParser(description= "compute sequence identity between two species")

arg_parser.add_argument("-b","--bedfile", help='.bed file in species 1 coordinates')
arg_parser.add_argument("-i","--iters", help='number of shuffles')
arg_parser.add_argument("-g","--genome_build", help='genome_build')
arg_parser.add_argument("-incl","--include", help='regions to include in shuffle')
arg_parser.add_argument("-o","--outdir", help='directory for files')

args = arg_parser.parse_args()

TEST_BED = args.bedfile
ITERS = args.iters
BUILD = args.genome_build
INCLUDE = args.include
SHUF_PATH = args.outdir

PATH, FILENAME, SAMPLE_ID = split_filename.split_filename(TEST_BED) 


###
#   functions
###


def loadConstants(build):  
    path_dict = {
                'hg19': ("/dors/capra_lab/users/fongsl/data/hg19_blacklist_gap_ensemblexon.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
                'hg38': ("/dors/capra_lab/users/fongsl/data/hg38_blacklist_gap_ensemblexon.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
                'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
                    }
    blacklist, sizes = path_dict[build]
    
    return blacklist, sizes


def shuffle(test_bed, shuf_path, sample_id, iter_, build, include):
    
    out = os.path.join(shuf_path, f"shuf-{sample_id}-{iter_}.bed")  # write file

    BLACKLIST, CHROM_SZ = loadConstants(build)

    if include is not None:    
        BEDshuf = f"bedtools shuffle -i {test_bed} -g {CHROM_SZ} -chrom -noOverlapping -maxTries 5000 -incl {include} -excl {test_bed} > {out}" 

    else:
        BEDshuf = f"bedtools shuffle -i {test_bed} -g {CHROM_SZ} -excl {BLACKLIST} -chrom -noOverlapping -maxTries 5000 > {out}" 

    #print(BEDshuf)
    os.system(BEDshuf)
    
    
###
#   Main
###


def main(argv):
    
    #num_cores = multiprocessing.cpu_count()
    num_cores = 16
    print("number of cores", num_cores)

    # run parallel jobs

    Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(shuffle)(TEST_BED, SHUF_PATH, SAMPLE_ID, i, BUILD, INCLUDE) for i in np.arange(int(ITERS)))


    shuffle(TEST_BED, SHUF_PATH, SAMPLE_ID, ITERS, BUILD, INCLUDE)
    
if __name__ == "__main__":
    main(sys.argv[1:])