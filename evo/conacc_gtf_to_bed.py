import argparse
import glob
import os
import pandas as pd
import subprocess
import sys

###
# args
###

arg_parser = argparse.ArgumentParser(description="post-run clean-up of conservation/acceleration scores w phyloP")

arg_parser.add_argument(
    "-msa", "--multiz", help='20-, 30-, 100-way multiz in hg38')

arg_parser.add_argument("-p", "--path", help = "path to dump results")


# PARSE THE ARGUMENTS
args = arg_parser.parse_args()

#IDX = args.index  # the index

MSA_WAY = args.multiz  # multiple sequence alignment.
DATA_PATH = args.path



query = os.path.join(DATA_PATH, "chr*", f"multiz{MSA_WAY}_br-*", f"chr*_conacc.bed")
FS = glob.glob(query)



def drop_gtf_cols(F):
    path = "/".join(F.split("/")[:-1])  # file name has full path. 
    
    os.chdir(path)  # go to the path
    
    cmd = f"cut -f 1,4,5,6,9 {F} > t && mv t {F}" # cut the files w/ info
    
    subprocess.call(cmd, shell = True)  # run cut command in shell. 

    
def format_id_col(F):
    
    f = pd.read_csv(F, sep = '\t', header = None)  # open file
    print(f.head())
    
    f[4] = f[4].apply(lambda x: x.split('"')[1])  # col 4 is the id col. String format id"bin_xxxxxx" to bin_xxxxxx
    
    f.to_csv(F, sep = '\t', header = False, index = False)  # save results 

    
def main(argv):    
    for F in FS:

        drop_gtf_cols(F)

        format_id_col(F)

if __name__ == "__main__":
    main(sys.argv[1:])