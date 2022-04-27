#!/usr/bin/env python
# coding: utf-8

# In[10]:


import configparser
import glob
import os
import subprocess
import sys



###
# args
###

arg_parser = argparse.ArgumentParser(description="post-run clean-up of conservation/acceleration scores w phyloP")

arg_parser.add_argument(
    "-chr", "--chromosome", help='20-, 30-, 100-way multiz in hg38')

arg_parser.add_argument("-p", "--path", help = "path to dump results")


# PARSE THE ARGUMENTS
args = arg_parser.parse_args()


CHR_PATH = args.chromosome
DATA_PATH = args.path



### count chr_regions

def concatnclean(path, chrnum, exp_count):
    
    # concat the output files
    outf = os.path.join(path, f"{chrnum}_conacc.bed")

    cat_cmd = f"cat {path}/{chrnum}_*_conacc.bed > {outf}"
    subprocess.call(cat_cmd, shell = True) 
    print("concat\n\n")
    
    num_lines_ = sum(1 for line in open(outf))

    if num_lines_ == exp_count:        
        cleanup_cmd = f"rm {path}/{chrnum}_*_conacc.bed"
        print("nclean\n\n")
        subprocess.call(cleanup_cmd, shell = True) 


# In[14]:



"""
# run a big for-loop 

per chromosome

get the original n lines of the bed file
count the lines of the resulting conacc files.
make sure that n lines of original chrN.bed == n lines resulting chrN_*_conacc.bed files

if n lines exp from chromosome file == n line result conacc file
concat chrN_*_conacc.bed files
clean up the chrN_*_conacc.bed files
"""


# In[16]:


def main(argv):
    query = os.path.join(CHR_PATH, "*.bed")

    chrs_ = glob.glob(query)

    for chr_ in chrs_:

        chrnum = (chr_.split("/")[-1]).split(".bed")[0]

        """
        # expected number of lines
        """

        exp_num_lines = sum(1 for line in open(chr_))

        # go count the result files, too 
        chr_dir = os.path.join(DATA_PATH, chrnum)

        query_res_paths = os.path.join(chr_dir, "multiz*")

        res_paths = glob.glob(query_res_paths) # get result paths
        print(res_paths)

        """
        # Obs number of lines per branch, neutral model
        """

        for res_path in res_paths:  # for every multiz results path

            outf = os.path.join(res_path, f"{chrnum}_conacc.bed")

            """
            check if already ran this cleanup in dir. 
            """
            if os.path.exists(outf) is True:

                print("\n\nalready run", chrnum, res_path)
                continue

            """
            # Count obs number of lines
            """
            query_res = os.path.join(res_path, f"{chrnum}_*_conacc.bed")  # get all the bed files

            res_ = glob.glob(query_res) 

            obs_num_lines = 0  # sum all the lines of the other files. 

            for res_f in res_:

                num_lines_ = sum(1 for line in open(res_f))  # sum all the lines of the other files. 

                obs_num_lines += num_lines_ 


            """
            # check if the results files are complete
            """


            if exp_num_lines == obs_num_lines:
                print("expected", chrnum, exp_num_lines, "obs", obs_num_lines)

                # complete! Clean up the directory. 
                concatnclean(res_path, chrnum, exp_num_lines)


            else:
                # not complete. Follow up
                print("\n\n not done yet", res_path, "expected",exp_num_lines, "obs",  obs_num_lines)

if __name__ == "__main__":
    main(sys.argv[1:])


# In[ ]:





