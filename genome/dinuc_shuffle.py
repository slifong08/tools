# sarahfong
# return dinuc scramble fasta from input fasta. 

import argparse
import numpy as np
from random import shuffle
import os, sys

def dinucScramble(seq, verbose=False):

    n=2 # pick size 2
    
    splits = [seq[i:i+n] for i in range(0, len(seq), n)]
    shuf_splits = splits[::]
    
    shuffle(shuf_splits)
    
    if verbose is True:
        print(np.unique(splits, return_counts=True))
        print(np.unique(shuf_splits, return_counts=True))
    

    return "".join(shuf_splits)

def main(argv):
    """scramble sequences in fasta"""

    ###
    # step 1 - parse input arguments
    ###
    
    arg_parser = argparse.ArgumentParser(description= "generate dinucleotide shuffles given list of sequences")

    arg_parser.add_argument("--input", help='input fasta or text file where each line is one sequence')
    args = arg_parser.parse_args()

    FA = args.input
    
    # file to write
    FA_SCRAMBLE = ".".join(FA.split(".")[:-1]) + ".scramble.fa"

    ###
    # step 2 - make dictionary of scrambled sequences from inputs
    ###
    
    # make fasta scramble dictionary
    fa_scrambles = {}

    # read and scramble fasta sequences
    with open(FA, "r") as fa_reader:
        seqid, counter=None, 0 # seq id, counter
        
        for row in fa_reader:
            row = row.strip("\n") 
            
            if ">" in row:  # handle the id
                seq_id = row + "_scramble"
                
            else:  # handle the sequence
                
                if seqid is None:  # incase this isn't a fasta file, but a list of sequences, 
                    seq_id=f">{counter}"
                    
                # scramble and store sequence in dictionary
                fa_scrambles[seq_id] = dinucScramble(row, verbose=False)  
                
                seqid=None # reset seqid
                
            counter +=1
    
    ###
    # step 3 - write scramble sequences
    ###
    
    with open(FA_SCRAMBLE, "w") as writer:
        
        for key, value in fa_scrambles.items():
            
            writer.write(f"{key}\n{value}\n")
            
        print("\n\nwrote", FA_SCRAMBLE)


if __name__ == "__main__":
    main(sys.argv[1:])