from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement, Seq
import numpy as np
import pandas as pd


def one_hot_encode(seq_vector, seq_length):
    """ One-hot encoding
        # handles N as 5th column
        from https://stackoverflow.com/questions/34263772/how-to-generate-one-hot-encoding-for-dna-sequences
    """
    # make vector of zeros
    X = np.zeros((seq_vector.shape[0], seq_length, 5))
    
    for idx in range(0, len(seq_vector)):
        
        seq = seq_vector[idx]
        mapping = dict(zip("ACGTN", range(5)))    
        seq2 = [mapping[i] for i in seq]
        try:
            X[idx] = np.eye(5)[seq2]
        except:
            print("error", idx, len(seq),  seq, seq2, len(seq2))

    return X

def fasta2dict(fasta_file):
    """Fasta 2 dictionary"""
    
    # Parse fasta
    fasta_dict = {}
    with open(fasta_file, "r") as reader:
        for value in SimpleFastaParser(reader):
            name, seq = value
            fasta_fwd[name] = seq

    return fasta_dict

def chrList():
    """return list of chromosomes"""
    
    chrs = []
    
    for n in np.arange(1,23):
        chrs.append(f"chr{n}")    
    # add sex chromosomes
    chrs.append("chrX")
    chrs.append("chrY")
    
    return chrs


def writeFa(seq_dict, out_fa):
    """write fa from sequence dictionary """
    
    with open(out_fa, "w") as writer:
        for seqid, seq in seq_dict.items():
            writer.write(f">{seqid}\n{seq}\n")


def splitChrTrainTestVal(chr_list, holdout_list=None):
    """return train (list), test (str), val (str) chromosomes
         if holdout_list is none, 
         randomly sample and hold out 1 chromosome for each test, val
     """

    # if heldout list not specified, 
    # randomly sample 2 held out chromosomes
    
    if holdout_list is None:

        holdout_list= list(np.random.choice(chr_list, 2))
    val_chr, test_chr = holdout_list
        
    # remove heldout chromosomes
    for i in holdout_list:
        chr_list.remove(i)

    return chr_list, val_chr, test_chr

def writeReverseComp(df):
    """return table w concatnetated reverse complement sequences, annotation column"""
    df["reverse"]=0
    dfr=df.copy()
    dfr["reverse"]=1 # re-annotate reverse complement
    dfr["seq"] = dfr["seq"].apply(lambda x: str(reverse_complement(Seq(x)))) # make reverse complement

    return pd.concat([df, dfr])
    