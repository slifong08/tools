#!/usr/bin/env python
# coding: utf-8


import argparse
import glob
import itertools
from joblib import Parallel, delayed
import multiprocessing
import os
import sys
import string
import subprocess
sys.path.append("/dors/capra_lab/users/fongsl/tools/py_/")
import split_filename

###
# args
###

arg_parser = argparse.ArgumentParser(description="estimate conservation/acceleration w phyloP")

#arg_parser.add_argument("-i", "index", help='array number', type=int)
arg_parser.add_argument("-chr", "--chromosome", help='e.g. chr12')
arg_parser.add_argument("-br", "--branches", help='hg38, rheMac8')
arg_parser.add_argument("-g", "--genome_build", help='hg38')
arg_parser.add_argument(
    "-msa", "--multiz", help='20-, 30-, 100-way multiz in hg38')
arg_parser.add_argument(
    "-mod", "--model", help='full", hg38-rheMac8', default="full")
arg_parser.add_argument("-o", "--outdirectory", help = "outdirectory to dump results")



# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

#IDX = args.index  # the index
CHRNUM = args.chromosome  # the chromosome (full path)
BRANCH = args.branches  # the branches to test.
MSA_WAY = args.multiz  # multiple sequence alignment.
MODEL = args.model
PHYLOP_DATA_PATH = args.outdirectory
GENOME_BUILD = args.genome_build


# CONSTANTS
RANDOM_SEED = 42
GROUP_SIZE = 10  # so that only 10 files are processed in parallel. 
# Keeps from overrunning memory requests on slurm


###
#  FUNCTIONS
###
def get_phylop_dict(msaway, genome_build):

    phylop_dict = {
        'phylop_bin':'/dors/capra_lab/bin/./phyloP', 
        'dors_maf_path': f'/dors/capra_lab/data/ucsc/{genome_build}',
        'maf':f'/dors/capra_lab/data/ucsc/hg38/multiz{msaway}/maf/{msaway}.maf.gz', 
        'branches':['hg38', 'rheMac8', 'hg38-rheMac8'], 
        'models' :['full', 'rheMac8_noOWM', 'hg38_noAPES'],
    }
    
    neutral_model_dict = {
                            'full': f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}.mod', 
                            'hg38-rhemac8': f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}_hg38-rheMac8.mod',
                            'rhemac8_noowm': f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}_rheMac8_noOWM.mod', 
                            'hg38_noapes':f'/dors/capra_lab/data/ucsc/{genome_build}/multiz{msaway}/{genome_build}.phastCons{msaway}_hg38_noAPES.mod'
                            }
    return phylop_dict, neutral_model_dict


# run phylop

def run_phylop(msaway, ocr, chrnum, path, random_seed, branch, model, model_dict, base, phylop):

    print(ocr, chrnum)
    n = ocr.split("-")[1]

    # neutral tree
    mod = model_dict[model]  # get the dictionary of the models

    # multiple sequence alignment file
    maf_zipped = os.path.join(base, f"multiz{msaway}", "maf", f"{chrnum}.maf.gz")
    maf_unzipped = maf_zipped.split(".gz")[0]

    # maf needs to be unzipped?
    if os.path.exists(maf_unzipped) is False:
        cmd = f"gunzip {maf_zipped}"
        subprocess.call(cmd, shell=True)

    # make outpath
    outpath = os.path.join(
        path, f"{chrnum}", f"multiz{msaway}_br-{branch}_mod-{model}")
    try:
        os.mkdir(outpath)
    except FileExistsError:
        pass

    # make outfile
    outf = os.path.join(outpath, f"{chrnum}_{n}_conacc.bed")

    # Already done phyloP analysis on this file?
    if os.path.exists(outf) is False or os.path.getsize(outf) == 0:

        # run phyloP!
        cmd = f"{phylop} --features {ocr} --msa-format MAF --method LRT --branch {branch} --mode CONACC -d {random_seed}         -g {mod} {maf_unzipped}> {outf}"
        print(cmd)
        # write run to log
        runlog_f = os.path.join(outpath, "runlog.txt")
        with open(runlog_f, "a") as runlog:
            runlog.write(cmd + "\n\n")

        # print(cmd)
        subprocess.call(cmd, shell=True)

        # check results
        if os.path.getsize(outf) > 0:

            # delete temp
            temp = os.path.join(path, f"temp_{chrnum}.bed")
            if os.path.exists(temp) is True:
                os.remove(temp)
                print("removed", temp)

        else:
            print("this didn't run", ocr)
    else:
        print("already processed", outf)

        # delete temp
        temp = os.path.join(path, f"temp_{chrnum}.bed")
        if os.path.exists(temp) is True:
            os.remove(temp)
            print("removed", temp)

    return outf

###
# # MAIN
###

def main(argv):
   
      
    PHYLOP_DICT, MODEL_DICT = get_phylop_dict(MSA_WAY, GENOME_BUILD)
   
    BASE = PHYLOP_DICT["dors_maf_path"]
    PHYLOP = PHYLOP_DICT["phylop_bin"] 
    CHR_DIR = os.path.join(PHYLOP_DATA_PATH, CHRNUM)
    os.chdir(CHR_DIR)  # change directory

    FS = os.listdir(CHR_DIR)  # get a list of the split files
    for F in FS:
        if "multiz" in F:
            FS.remove(F)
        elif "stderr" in F:
            FS.remove(F)

    # prepare to run parallel jobs as

    num_cores = GROUP_SIZE
    print("number of cores", num_cores, multiprocessing.cpu_count())

    # run parallel jobs

    Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(run_phylop)(
        MSA_WAY, ocr, CHRNUM, PHYLOP_DATA_PATH, RANDOM_SEED, BRANCH, MODEL, MODEL_DICT, BASE, PHYLOP) for ocr in FS)

if __name__ == "__main__":
    main(sys.argv[1:])