import argparse
import glob
import os, sys

"""
GOAL: Transform .bed file into fasta, chr by chr using bedtools "getfasta" command. 

NOTES:
- human DNA only
- writes .fa file to same directory as .bed file


"""

###
#   arguments
###

arg_parser = argparse.ArgumentParser(description="turn .bed into .fa, chromosome by chromosome")

arg_parser.add_argument("bed", help='.bed file w/ full path')
arg_parser.add_argument("-b", "--build", type=str, default='hg19', choices=['hg19', 'hg38', "hs1"],
                        help='human genome assembly; default=hg19')

args = arg_parser.parse_args()

TEST_BED = args.bed
BUILD = args.build


### 
# Functions
### 

def get_path_filename_sample_id(test_bed):
    
    test_path, filename = os.path.split(test_bed)
    
    sample_id = os.path.splitext(filename)[0]

    return test_path, filename, sample_id

    
    
def get_dna(genome_build):
    
    if "capra" in os.getcwd():  # accre
        dna_path = f"/dors/capra_lab/data/dna/human/{genome_build}"
    else:  # wynton
        dna_path =f"/wynton/home/ahituv/fongsl/dna/{genome_build}"
        
    return dna_path



###
#   main
###

def main(argv):

    # turn bed into fasta using bedtools getfasta command

    test_path, filename, sample_id = get_path_filename_sample_id(TEST_BED)
 
    # path to input DNA fasta files 
    dna_path = get_dna(BUILD)
    
    # file to read, write
    fasta_in = os.path.join(dna_path, f"{BUILD}.fa")
    fasta_out = os.path.join(test_path, f"{sample_id}.fa")
    
    # unzip fasta_in if zipped. 
    if os.path.exists(f"{BUILD}.fa.gz") is True:
        os.system(f"gunzip {BUILD}.fa.gz")
    
    # intersect test.bed w/ fasta. 
    cmd = f"bedtools getfasta -fi {fasta_in} -bed {TEST_BED} -fo {fasta_out}"

    os.system(cmd)
    
    print("wrote", fasta_out)
    
    
if __name__ == "__main__":
    main(sys.argv[1:])

