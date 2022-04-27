import argparse
from collections import Counter
import os

"""
Output: .tsv with sequence id,  GC count, and GC density

Input: .fa file

Notes
- I recommend running fasta_from_bed.py to make fasta files for human sequences. 

"""

###
#   arguments
###

arg_parser = argparse.ArgumentParser(description="quantify GC content (count nucleotides, compute density")

arg_parser.add_argument("f","--fasta", help='.fa file w/ full path')

args = arg_parser.parse_args()

FASTA = args.fasta


### 
# Functions
### 

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
        
        
def format_fa(fasta_file):
    
    path, filename= os.path.split(fasta_file)
    split = os.path.splitext(filename)[0]

    outf = os.path.join(path, f"{split}_GC.tsv" ) # path to write
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

    out_fasta = format_fa(FASTA)
    
    
if __name__ == "__main__":
    main(sys.argv[1:])


        

