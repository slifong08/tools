import glob
import numpy as np
import os, sys
import subprocess


###
# functions
###


def split_into_chr_bed(data, path):

    os.chdir(path)
    n_files = len(glob.glob("*.bed"))
    if n_files!=22:
        cmd = ''' awk '{print>$1".bed"}' %s''' %  data
        subprocess.call(cmd, shell = True)
        print(cmd)
        
        cmd = "rm *random*.bed"
        subprocess.call(cmd, shell = True)

        cmd = "rm \#chr.bed"
        subprocess.call(cmd, shell = True)
        
        
    files = glob.glob("*.bed")

    for n in files:
        if os.path.getsize(n)> 0:
            if "regions" in data:
                cols = '1,2,3,6'
            else:
                cols = '1,2,3,4'
                
            cmd = f'cut -f {cols} {n} > t'
            subprocess.call(cmd, shell = True)  ## keep coordinates, bin annotation
            
            cmd = f'mv t {n}'
            subprocess.call(cmd, shell = True)  # rename the file. 

def make_chr_list():
    n = list(np.arange(1, 23))
    #n.append("X")

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list