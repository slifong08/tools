import argparse
from collections import Counter
import glob
import os, sys

"""
Output: .tsv with sequence id,  GC count, and GC density

Input: .bed file

Notes
- turns .bed -> .fa, then evaluates GC
- human DNA only
- writes .fa file to same directory as .bed file
"""

###
#   arguments
###

arg_parser = argparse.ArgumentParser(description="turn .bed into .fa, chromosome by chromosome, get GC content")

arg_parser.add_argument("-i", "--bed", help='.bed file w/ full path')
arg_parser.add_argument("-b", "--build", type=str, default='hg19', choices=['hg19', 'hg38'],
                        help='human genome assembly; default=hg19')
arg_parser.add_argument("-o", "--outdirectory", help='directory to dump gc information')




args = arg_parser.parse_args()

TEST_BED = args.bed
BUILD = args.build
OUTDIR = args.outdirectory

### 
# Functions
### 

def get_path_filename_sample_id(test_bed):
    
    test_path, filename = os.path.split(test_bed)
    
    sample_id = os.path.splitext(filename)[0]

    return test_path, filename, sample_id
    
    
def split_by_chr(test_bed, outdir):
        
        test_path, filename, sample_id = get_path_filename_sample_id(test_bed)
        
        # go to bed dir
        os.chdir(outdir)
        
        # split test into chrN.bed files
        chr_cmd = "awk '{print >$1\"_%s_temp.bed\"}' %s" % (sample_id, test_bed) 
        print(chr_cmd)
        os.system(chr_cmd)
    
        # glob chr split bed files. 
        bed_chr_list = glob.glob(os.path.join(outdir, f"chr*_{sample_id}_temp.bed"))  # glob chromosomes

        return bed_chr_list
    
    
def get_dna(genome_build):
    
    dna_path = f"/dors/capra_lab/data/dna/human/{genome_build}"
    
    return dna_path


def getfasta(bed_chr_list, test_bed, genome_build, outdir):

    test_path, filename, sample_id = get_path_filename_sample_id(test_bed)
    
    # file to write
    out_fasta = os.path.join(outdir, f"{sample_id}.fa")


    # intersect test chrN.bed w/ fasta. 
    for bed_chr in bed_chr_list: 

        chr_num = os.path.basename(bed_chr).split("_")[0]

        # path to input DNA fasta files 
        dna_path = get_dna(genome_build)

        fasta_in = os.path.join(dna_path, f"{genome_build}.fasta")
        fasta_out = os.path.join(outdir, f"{sample_id}_{chr_num}.fa")

        cmd = f"bedtools getfasta -fi {fasta_in} -bed {bed_chr} -fo {fasta_out}"
        #print(cmd)
        os.system(cmd)

        # clean up chr temp file
        os.remove(bed_chr)

    # concat all the chromosomes
    cat_handle = os.path.join(outdir, f"{sample_id}_chr*.fa")
    cat_cmd = f"cat {cat_handle} > {out_fasta}"
    os.system(cat_cmd)

    # remove chromosome fasta files. 
    cleanup_cmd = f"rm {cat_handle}"
    os.system(cleanup_cmd)

    return out_fasta

### 
# GC Functions
### 

def countgc(sequence):
    gc = []
    total = []
    letters = ["G", "C", "g", "c"]
    counts = Counter(sequence)

    for letter in letters:
        gc.append(int(counts[letter]))  # count all the Gs and all the Cs

    gc_sum = sum(gc)
    gc_density = gc_sum/len(sequence) 
    
    return gc_sum, gc_density


def writerow(outf, rows):
     with open(outf, "w") as writer:
        for row in rows:
            writer.write(row)
        writer.close()
        
        
def format_fa(fasta_file, outdir):
    
    path, filename= os.path.split(fasta_file)
    split = os.path.splitext(filename)[0]

    outf = os.path.join(outdir, f"{split}_GC.tsv" ) # path to write
    print(outf)
    already_processed = []
    
    with open(fasta_file, "r") as ff:
    
        fasta = ff.readlines()

        new_line = []  # for writing new line
        rows = []
        
        for n, line in enumerate(fasta):
            #print(line)
            
            if ">" in line:  # this is a new sequence
                seq_id = (line.split(">")[1]).strip("\n")
                
            elif line == "":
                seq_id = None
                continue
            else:  # this is the sequence
                sequence = line
                
                gc_count, gc_density = countgc(sequence) # count GC content, density for sequence

                row = f"{seq_id}\t{gc_count}\t{gc_density}\n"
                
                if seq_id not in already_processed:
                    
                    # append seq id to already_processed list
                    already_processed.append(seq_id) 
                    
                    # append row
                    rows.append(row)
                else:
                    print(f"\nCHECK .FA REDUNDANCY. Seq_id {seq_id} is already in dictionary.\n")
                    

        writerow(outf, rows)
                   
    return outf

def main(argv):
                          
    # split bed file by chromosome
    # turn bed into fasta using bedtools getfasta command
    bed_chr_list = split_by_chr(TEST_BED, OUTDIR)
    
    FASTA = getfasta(bed_chr_list, TEST_BED, BUILD, OUTDIR)
   
    print(FASTA)
    
    out_fasta = format_fa(FASTA, OUTDIR)
    
    
if __name__ == "__main__":
    main(sys.argv[1:])
        