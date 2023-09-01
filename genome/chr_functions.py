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

def makeFullChrList():
    n = list(np.arange(1, 23))
    n.extend(["X", "Y"])

    chr_list = []
    for num in n:
        chrn = "chr" + str(num)
        chr_list.append(chrn)

    return chr_list

def makeCoorAnnot(df, chr_colname, start_colname, end_colname, id_name):
    """
    make a bed coordinates column in dataframe
    
    input 
        df (pd dataframe) - dataframe to transform
        chr_colname (str) - name of #chr column
        start_colname (str) - name of start column
        end_colname (str) - name of end column
        id_name (str) - name to call new coordinate column

    method
        make coordinate column from #chr, start, end columns

    return
        df (pd dataframe) - dataframe with one extra column for the coordinate annotation
        
        
    """
    df[f"{id_name}.coor"] = df[chr_colname].map(str) +  ":"  + df[start_colname].map(str) + "-" + df[end_colname].map(str)
    
    return df