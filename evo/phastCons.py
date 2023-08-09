#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import subprocess
import sys

sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
sys.path.append("/dors/capra_lab/users/fongsl/tools/genome/")

import config_readwrite as crw
import chr_functions
import split_filename

"""
Intersect .bed file w/ phastcons elements 

Output: .bed file w/ phastcons overlap count (i.e. how many PhastCons elements overlap. 

Input: 
    - .bed: file, 
    - genome build: str
    - msa-way: str, must include 'way' (e.g. '30way')


Notes


"""

###
#   arguments
###

arg_parser = argparse.ArgumentParser(description= "compute sequence identity between two species")

arg_parser.add_argument("-b","--bedfile", help='.bed file in species 1 coordinates')
arg_parser.add_argument("-g","--genome_build", help='e.g. hg38')
arg_parser.add_argument("-m","--msa", help='multiz: 100way, 30way, 20way')
arg_parser.add_argument("-o","--outdirectory", help='out directory to save results')


args = arg_parser.parse_args()

TEST_BED = args.bedfile
GENOME_BUILD = args.genome_build
MSAWAY = args.msa
OUTDIR = args.outdirectory


# # function: intersect test with phastCons 

def load_phastcons(msaway, genome_build):

    return f"/dors/capra_lab/data/evolutionary_conservation/phastcons/{genome_build}/phastConsElements{msaway}_{genome_build}.bed"


def bed_intersect(test_bed, phast, msaway, outdir):

    path, filename, sample_id =  split_filename.split_filename(test_bed) 
    
    outfile = os.path.join(outdir, f"{sample_id}_phastcons-{msaway}.bed")  # make out file

    """
    # run bedtools command, 
    # -c count the number of phastcons overlaps, which
    # reports 0 if no overlaps w/ Phastcons

    """
    cmd = f"bedtools intersect -a {test_bed} -b {phast} -c > {outfile}"  

    print("intersecting", test_bed, "with", phast)
    
    subprocess.call(cmd, shell = True)  # run it. 

    
    return outfile

# # run 

# In[ ]:

def main(argv):
    
    phast = load_phastcons(MSAWAY, GENOME_BUILD)
    
    bed_intersect(TEST_BED, phast, MSAWAY, OUTDIR)

if __name__ == "__main__":
    main(sys.argv[1:])






