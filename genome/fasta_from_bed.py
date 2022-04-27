import argparse
import glob
import os

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
arg_parser.add_argument("-b", "--build", type=str, default='hg19', choices=['hg19', 'hg38'],
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
    
    
def split_by_chr(test_bed):
        
        test_path, filename, sample_id = get_path_filename_sample_id(test_bed)
        
        # go to bed dir
        os.chdir(test_path)
        
        # split test into chrN.bed files
        chr_cmd = "awk '{print >$1\"_%s_temp.bed\"}' %s" % (sample_id, test_bed) 
        os.system(chr_cmd)
    
        # glob chr split bed files. 
        bed_chr_list = glob.glob(os.path.join(test_path, f"chr*_{sample_id}_temp.bed"))  # glob chromosomes
        
        return bed_chr_list
    
    
def get_dna(genome_build):
    
    dna_path = f"/dors/capra_lab/data/dna/human/{genome_build}"
    
    return dna_path


def getfasta(bed_chr_list, test_bed, genome_build):

    test_path, filename, sample_id = get_path_filename_sample_id(test_bed)
    
    # file to write
    out_fasta = os.path.join(test_path, f"{sample_id}.fa")



    # intersect test chrN.bed w/ fasta. 
    for bed_chr in bed_chr_list: 

        chr_num = os.path.basename(bed_chr).split("_")[0]

        # path to input DNA fasta files 
        dna_path = get_dna(genome_build)

        fasta_in = os.path.join(dna_path, f"{genome_build}.fasta")
        fasta_out = os.path.join(test_path, f"{sample_id}_{chr_num}.fa")

        cmd = f"bedtools getfasta -fi {fasta_in} -bed {bed_chr} -fo {fasta_out}"
        #print(cmd)
        os.system(cmd)

        # clean up chr temp file
        os.remove(bed_chr)

    # concat all the chromosomes
    cat_handle = os.path.join(test_path, f"{sample_id}_chr*.fa")
    cat_cmd = f"cat {cat_handle} > {out_fasta}"
    os.system(cat_cmd)

    # remove chromosome fasta files. 
    cleanup_cmd = f"rm {cat_handle}"
    os.system(cleanup_cmd)

    return out_fasta


###
#   main
###

def main(argv):

    # split bed file by chromosome
    # turn bed into fasta using bedtools getfasta command
    bed_chr_list = split_by_chr(TEST_BED)
    out_fa = getfasta(bed_chr_list, TEST_BED, BUILD)
    
    
if __name__ == "__main__":
    main(sys.argv[1:])

