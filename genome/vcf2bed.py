"""
sarahfong

convert 1-index .vcf to 0-index .bed 

unzipped .vcf only

"""

import os, sys

def vcf2bed(vcf):
    """
    turn .vcf -> bed
    
    input
        vcf (str) - vcf file name
        
    method
        1. make bed output file
        2. open vcf to read, bed to write
        3, write column names
        4. parse through vcf lines
        5. get VCF info, convert 1-based position to 0-based start
        6. write new line as tab separated string. 
        
    return 
        bed (str) - output bed file name
    """
    
    #1
    BED = vcf.strip(".vcf") + ".bed"
    
    #2
    reader, writer = open(vcf, "r"), open(BED, "w")

    #2 write column names
    col_names =["#CHROM","POS-1","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    writer.write("\t".join(col_names)+"\n")
    
    #4
    for line in reader:
        if "#" not in line:  # skip info lines
            
            #5
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.strip("\n").split("\t")
            START = int(POS)-1  # 1-index converts to zero-index
            
            #6
            new_line = f"chr{CHROM}\t{START}\t{POS}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\n"
            writer.write(new_line)
        
    writer.close()
    
    return BED