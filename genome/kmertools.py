"""
tools I wrote to handle kmers. 
"""
import numpy as np

def CountKmers(sequence, windowsize, stepsize):

    """
    Count unique kmers from sequence; break up sequence into equally sized kmers w sliding window of with windowsize, stepsize
    
    required packages
        - numpy
        
    inputs
        sequence (str) - sequence to break into kmers
        windowsize (int) - windowsize to make for sequence
        stepsize (int) - steps between windows
        
    method
        1. calculate total n possible kmers given sequence length and window size. 
        2. instantiate collection list for kmers
        3. iterate through sequence w/ stepsize, starting at index=zero
        4. if index is less than or equal to n possible kmers, append sequence window to kmer list
        5. count the numbers of unique kmers and obeserved counts 
    
    return 
        vals (np.array) - array of unique kmers
        counts (np.array) - array of counts corresponding to unique kmers in values. 
        
        
        
    """
    #1 n possible kmers
    npossible = (len(sequence)-windowsize)

    #2
    kmers=[]
    
    #3
    for n in range(0, len(sequence),stepsize):
        
        #4
        if n<=npossible: # append sequence windows within range of possible windows
        
            kmers.append(sequence[n:n+windowsize])  # append sequence window to kmer list
            print(sequence[n:n+windowsize], n)
    
    #5
    vals, counts = np.unique(kmers, return_counts=True) # count window-size nucleotide content. 
    
    return vals, counts