from Bio.SeqUtils import gc_fraction, GC

def count_dinucleotide(sequence):
    
    """
    input 
        sequence (str)
    method
        1. instantiate GC count variable, last_dinucleotide variable. Upper() sequence
        2. iterate, step-size=1, through sequence, counting GC/CG
        3. if last dinucleotide pair was not in the previous step, count this dinucleotide. 
            Set the last_dinuc step at this index. 
        4. if the last dinucleotide pair WAS adjacent to this step, do not count, continue. 
        5. estimate fraction dinucleotide 
    
    return count, fraction of GC dinucleotide
        
    """
    #1
    gc_dinuc_count, last_dinuc = 0, 0
    sequence = sequence.upper()
    
    #2
    for n in range(len(sequence)+1):  

        seq=sequence[n:n+2]  #get index of sequence and next nucleotide
        
        if seq=="GC" or seq=="CG":
            
            # 3
            if last_dinuc!=n-1:  # if last dinucleotide was not the previous step
                gc_dinuc_count +=1
                last_dinuc = n
                #print("Count sequence dinucleotide")

            # 4  handle cases where gc isn't dinucleotide, but trinucleotide
            elif last_dinuc ==n-1: 
                #print("last sequence was dinucleotide")
                continue
    #5
    gc_dinuc_frac = gc_dinuc_count/len(sequence)
    
    return gc_dinuc_count, gc_dinuc_frac


def count_gc(sequence):
    """
    count the frequency of G's and C's, the dinucleotide frequency
    
    input
        
        seq (str) - sequence to count GC frequency
    
    method
        1. quantify %GC with SeqUtils.gc_fraction function
        
    return 
        gc_fraction (float) - %gc
        
    """
 
    return GC(sequence), gc_fraction(sequence)