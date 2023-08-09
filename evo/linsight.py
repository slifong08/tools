"""
### Summary ###


# intersect .bed with linsight per chromosome.
# (linsight data is >1gb, so I split by chromosome to handle smaller data)
# ultimately, output file adds 5 columns to original .bed file
# 5 columns are ["chr_lin",  "start_lin", "end_lin","linsight_score", "overlap"]


NOTES 
- hg19 only!!
"""

# load modules

import argparse
import glob
import numpy as np
import os, sys
import subprocess


# PATHS AND FILES

arg_parser = argparse.ArgumentParser(description="Calculate linsight for bedfile.")

arg_parser.add_argument("bedfile", help='bed file w/ full path')

args = arg_parser.parse_args()

ENHANCERF = args.bedfile
ENHANCERPATH, ENHANCERFILE = os.path.split(ENHANCERF)
ENHANCER_ID = ".".join(ENHANCERFILE.split(".")[:-1])


LINSIGHTPATH = "/dors/capra_lab/data/evolutionary_conservation/linsight/"
LINSIGHTFILE = "LINSIGHT.bed.gz"
LINSIGHTF = os.path.join(LINSIGHTPATH, LINSIGHTFILE)
LINSIGHT_ID = "linsight"

OUTPATH = os.path.join(ENHANCERPATH, "linsight")

if not os.path.exists(OUTPATH):
    os.makedirs(OUTPATH)
    print("Directory '%s' created" %OUTPATH)


#%% FUNCTIONS


def chr_lst(): # make a list of autosomes

    chr_list = []

    for n in np.arange(1,23):

         chr_list.append("chr"+str(n))

    return chr_list


#split up the linsight/enhancer file by chromosome.
def split_by_chr(path, f, id_):

    chr_fs = glob.glob(os.path.join(path, f"chr*_{id_}.bed")) # check that you haven't split already.

    if len(chr_fs) == 0:

        os.chdir(path)

        split_cmd = "awk '{print >$1\"_%s.bed\"}' %s" % (id_, f)

        subprocess.call(split_cmd, shell = True)
        print(split_cmd, "\n", f, id_)

        print("split", f)

    else:

        print("already split", f)


def make_chr_file_dict(chr_list, path, id_): # make a dictionary of chr.

    file_dict = {}

    for chr_num in chr_list:

        f_chr = os.path.join(path, f"{chr_num}_{id_}.bed")

        if id_ == "linsight": # linsight files are zipped.

            f_chr = f_chr + ".gz"

        file_dict[chr_num] = f_chr

    return file_dict


def bed_intersect(enh_chr, lin_chr, outfile):

        bed_cmd = f"bedtools intersect -a {enh_chr} -b {lin_chr} -wao > {outfile}"
        print(bed_cmd)
        subprocess.call(bed_cmd, shell = True)
        
        print("finished", outfile)


def cleanup_dir(path, id_, outpath):

        # concat chromosome files
        CATF = os.path.join(outpath, f"{id_}_x_linsight.bed")
        cat_query = os.path.join(outpath, "chr*_x_linsight.bed")
                             
        cat_cmd = f"cat {cat_query} > {CATF}"
        subprocess.call(cat_cmd, shell = True)
        print(cat_cmd)
                             
        # clean up the chromosome files.
        cleanup_cmd = f"rm {path}chr*_{id_}.bed"
        subprocess.call(cleanup_cmd, shell = True)

        cleanup_cmd = f"rm {outpath}/chr*.bed"
        subprocess.call(cleanup_cmd, shell = True)


        return CATF


def expand_per_basepair_score(df): # expand linsight score per basepair.

    expanded_scores = np.repeat(df.linsight_score, df.lin_len)

    return expanded_scores


#%% run functions


chr_list = chr_lst() # generate list of chromosomes

# make linsight chr files, dictionary. I did this in 2019.
split_by_chr(LINSIGHTPATH, LINSIGHTFILE, LINSIGHT_ID)
linchr_dict = make_chr_file_dict(chr_list, LINSIGHTPATH, LINSIGHT_ID)

#%% make enhancer chr files, dictionary
split_by_chr(ENHANCERPATH, ENHANCERFILE, ENHANCER_ID)
enhchr_dict = make_chr_file_dict(chr_list, ENHANCERPATH, ENHANCER_ID)


#%%


# bedtools intersection of enhancer x linsight, writing all the overlaps.
for chr_num in chr_list:

    enhancer_chr = enhchr_dict[chr_num]
    linsight_chr = linchr_dict[chr_num]

    OUTFILE = chr_num + "_x_linsight.bed"
    OUTF = os.path.join(OUTPATH, OUTFILE) # write outfile

    bed_intersect(enhancer_chr, linsight_chr, OUTF)


#%% clean up the mess

intersection_file = cleanup_dir(ENHANCERPATH, ENHANCER_ID, OUTPATH)
