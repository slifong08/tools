from random import shuffle
import numpy as np

def dinucScramble(seq, n=1, verbose=False):
    """ return one sequence that has been dinucleotide shuffled (when n=1) | trinucleotide shuffle (when n=2). 
        if verbose is True, will report the original and shuffled di/trinucleotide counts. These should be identical.  
    """
    
    splits = [seq[i:i+n] for i in range(0, len(seq), n)] # make a list of nucleotide pairs
    shuf_splits = splits[::] # copy the nucleotide pair list
    
    shuffle(shuf_splits)  # shuffle the nucleotide pairs
    
    if verbose is True:
        print(np.unique(splits, return_counts=True))
        print(np.unique(shuf_splits, return_counts=True))
    

    return "".join(shuf_splits)  # return a joined 

