from Bio.SeqIO.FastaIO import SimpleFastaParser
import os, sys
from pathlib import Path

FASTA = sys.argv[1]

#%% FUNCTIONS

def readFastaWriteBed(fasta, outbed):
    
    # touch the file
    if os.path.exists(outbed) is False:
        Path(outbed).touch()

    #parse fasta file
    with open(fasta, "r") as handle:
        for values in SimpleFastaParser(handle):
            
            coor, seq = values
            
            # if coordinates are annotated, str-split to get information
            if "chr" in coor and "-" in coor:
                # get bed coordinates
                chr_ = coor.split(":")[0]
                start = (coor.split(":")[1]).split("-")[0]
                end = (coor.split(":")[1]).split("-")[1]
                
            
                # get line of bedfile 
                bedline = f"{chr_}\t{start}\t{end}\n"
                #print(bedline)
            
                # concat command
                os.system(f"echo {bedline} >> {outbed}")

            else:
                continue
            
    handle.close()
        
#%%

def main(argv):

    outbed = os.path.splitext(FASTA)[0] + ".bed"

    readFastaWriteBed(FASTA, outbed)

    
if __name__ == "__main__":
    main(sys.argv[1:])