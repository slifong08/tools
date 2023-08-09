"""
From bed file w DNA sequence str, calculate GC count, GC %, GC dinucleotide content

Input: 
    .bed file (str) - full path to bed file
    -b (str) - genome build. Default is hg38
    -o (str) - full path to output directory

Method 
    1. instantiate and parse arguments
    2.  get file names, path, sample id
    3. get chr fa files
    4. split bedfile into chromosomes
    5. turn bed into fasta using bedtools getfasta command
    6. combine all chr fastas together
    7. get GC fraction per bed file
    8. move fasta files to an fa directory


Output: 
    .tsv with sequence id,  GC count, and GC density

Notes
- turns .bed -> .fa, then evaluates GC
- human DNA only
- writes .fa file to same directory as .bed file
"""

import argparse

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import gc_fraction

import gc_content as gc
import glob
from joblib import Parallel, delayed
import os, sys

import subprocess as sp



###
#   arguments
###

arg_parser = argparse.ArgumentParser(description="turn .bed into .fa, chromosome by chromosome, get GC content")

arg_parser.add_argument("bed", type=str, help='.bed file w/ full path')
arg_parser.add_argument("-b", "--build", type=str, default='hg38', choices=['hg19', 'hg38'],
                        help='human genome assembly; default=hg38')
arg_parser.add_argument("-o", "--outdirectory", type=str, help='directory to write gc information')

# parse the arguments 
args = arg_parser.parse_args()

TEST_BED = args.bed
BUILD = args.build
OUTDIR = args.outdirectory

### 
# Functions
### 

def splitFile(test_file):
    """
    split test file into path, filename, and sample_id
    
    requirements
        os 
        
    input
        test_file (str) - full path + test_file name

    method
        1. path, file = os.splitpath
        2. filename = os.splittext
    
    return 
        path, file, filename (strs)
    """
    
    path, file = os.path.split(test_file)
    
    filename = os.path.splitext(file)[0]

    return path, file, filename
    
    
def splitByChr(test_bed, filename, outdir):
        """
        Split a file by chromosome name in awk 
        
        input
            test_bed (str) - .bed file w/ full path to split by chromosome
            outdir (str) - path to write output files
            
        method 
    
            1. change directory to outdirectory
            2. awk split on chromosome column
            3. glob resultant files.
        
        return 
            bed_chr_list (list) - list of chr-split .bed files
        
        """
        
        
        if os.getcwd()!=outdir:
            #1 go to bed dir
            os.chdir(outdir)
        
        #2 split test into chrN.bed files
        
    
        # old awk
        cmd = "awk '{print >$1\"_%s.bed\"}' %s" % (filename, test_bed) 
        print(cmd)
        os.system(cmd)
        #args = ["awk", '{print >$1\"_%s_temp.bed\"}' % filename, test_bed]
        #p = sp.Popen(args, stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE)
        #print(p.stdout.readline()) # read the first line

        #3 glob chr split bed files. 
        bed_chr_list = glob.glob(os.path.join(outdir, f"chr*_{filename}.bed"))  # glob chromosomes

        return bed_chr_list
    
    
def getDNA(genome_build):
    
    """
    return path to chr.fa.gz files
    """
    
    if "capra" in os.getcwd():    
        dna_path = f"/dors/capra_lab/data/dna/human/{genome_build}"

    else:
        dna_path = f"/wynton/home/ahituv/fongsl/dna/{genome_build}/chromosomes"
    
    return dna_path


def getFasta(bed_chr, dna_path, outdir):
    """
    get fasta from bed per chromosome 
    
    requirement
        bedtools getfasta
    
    input
        bed_chr (str) - chr.bed file w full path
        dna_path (str) - path to chr.fa files
        outdir (str) - path to write .fa outfile to

    method
        1. str split for filename, chr_num, 
        2. instantiate fasta input and output files, zipped input file
        3. unzip source chr.fa file if zipped. 
        4. use bedtools getfasta to get fasta from .bed
        5. clean up the bed file
        
    """
    #1
    filename = (bed_chr.split("/")[-1]).split(".bed")[0]
    chr_num = filename.split("_")[0]
   
    #2
    fasta_in = os.path.join(dna_path, f"{chr_num}.fa") 
    fasta_out = os.path.join(outdir, f"{filename}.fa")  # write
    
    zipped_in = fasta_in + ".gz"
   
    #3
    if os.path.exists(zipped_in) is True:
        print("unzipping", zipped_in)
        sp.call(f"gunzip {zipped_in}", shell=True)

    
    #4 
    cmd = f"bedtools getfasta -fi {fasta_in} -bed {bed_chr} -fo {fasta_out}"
    
    if os.path.exists(fasta_out) is False:
        sp.call(cmd, shell=True)
    
        line_count = sum(1 for line in open(fasta_out, "r").readlines())

        #5 clean up chr temp file
        if line_count > 0:
            os.remove(bed_chr) 

            
def combineChrs(outdir, filename):
    """
    cat chr.fa files and delete individual files
    """

    # file to write
    out_fasta = os.path.join(outdir, f"{filename}.fa")

    # concat all the chromosomes
    cat_handle = os.path.join(outdir, f"chr*{filename}*.fa")
    
    os.system(f"cat {cat_handle} > {out_fasta}")

    # remove chromosome fasta files. 
    os.system(f"rm {cat_handle}")
    
    
    return out_fasta

### 
# GC Functions
### 


def countGC(sequence):
    """
    return gc_fraction, gc_dinucleotide content fraction (custom script) 
    """

    # dinucleotide info
    dinuc_count, dinuc_frac = gc.count_dinucleotide(sequence)
    
    return gc_fraction(sequence), dinuc_frac


def writeRow(outf, rows):
    
    """
    write to outfile list of rows, one row per index. 
    """
    
    with open(outf, "w") as writer:
        for row in rows:
            writer.write(row)
        writer.close()
        
def gcFromFasta(fasta_file, filename, outdir):
    """
    format the fasta file into GC content
    
    require 
        Bio.SeqIO.FastaIO.SimpleFastaParser
        
    input
        fasta_file (str) - path to fasta 
        filename (str)-name of file
        outdir (str) - path to write to 
        
    method 
        1. make output file
        2. instantiate rows, already_processed lists
        3. parse through fasta w/ SimpleFastaParser
            get seqid, sequence
        4. calculate gc content (fraction G|C, fraction G+C)
        5. append data to lists
        6. write gc content information to file

    return 
        outfile (str) - written gc content file. 
    """
    
    #1
    outf = os.path.join(outdir, f"{filename}_GC.tsv" ) # path to write
    print("writing GC content", outf)
    
    #2
    rows, already_processed = [], []  # collect the sequence ids that have already been processed. 
 
    #3
    with open(fasta_file, "r") as ff:
        for values in SimpleFastaParser(ff):
            seq_id, sequence = values 
            
            #4
            gc_frac, gc_dinuc_frac = countGC(sequence.upper()) # count GC content, density for sequence

            row = f"{seq_id}\t{gc_frac}\t{gc_dinuc_frac}\n"
                
            if seq_id not in already_processed:
                    
                #5 append row to rows, seq id to already_processed lists
                rows.append(row), already_processed.append(seq_id) 
                        
            else:
                print(f"\nCHECK .FA REDUNDANCY. Seq_id {seq_id} is already written.\n")
                    
    #6
    writeRow(outf, rows)
                   
    return outf
    
def moveFa(outdir, fasta_file):
    """
    move fasta_file to a separate fa directory
    """
    fa_dir = os.path.join(os.path.dirname(outdir), "fa") 
    
    if os.path.exists(fa_dir) is False:  # make sure the fa directory exists first
        os.mkdir(fa_dir)
    
    # move the input fasta file to the fa directory
    cmd = f"mv {fasta_file} {fa_dir}"
    
    os.system(cmd)
    
    
def main(argv):
    

    #2 str split file name
    path, file, filename = splitFile(TEST_BED)
    
    #3 path to input DNA fasta files 
    dna_path = getDNA(BUILD)
    
    #4 split bed file by chromosome
    bed_chr_list = splitByChr(TEST_BED, filename, OUTDIR)

    # run get fasta per chromosome. 
    num_cores = len(bed_chr_list)

    #5 turn bed into fasta using bedtools getfasta command
    Parallel(n_jobs=24, verbose=100, prefer="threads")(delayed(getFasta)(i, dna_path, OUTDIR) for i in bed_chr_list)

    #6 combine all chr fastas together
    FASTA = combineChrs(OUTDIR, filename)

    #7 get GC fraction per bed file
    out_fasta = gcFromFasta(FASTA, filename, OUTDIR)
    
    #8 move fasta files to an fa directory
    moveFa(OUTDIR, FASTA)

if __name__ == "__main__":
    main(sys.argv[1:])