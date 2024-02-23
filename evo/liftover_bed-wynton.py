import argparse
import os, sys
import subprocess


arg_parser = argparse.ArgumentParser()


arg_parser.add_argument("bedfile",  
                        type=str, help='bed file w/ full path'
                       )
arg_parser.add_argument("from_build", 
                        type=str, help='genome build of the original file')

arg_parser.add_argument("to_build", 
                        type=str, help='genome build to liftover to')

arg_parser.add_argument("-c", "--chain", 
                        required=False, 
                        type=str, help='genome build to liftover to')
arg_parser.add_argument("-m", "--minmatch", 
                        default="0.95",
                        type=str, help='min match param from liftOver, \
                        "Minimum ratio of bases that must remap", default = 0.95')

args = arg_parser.parse_args()

# args
F, FROM, TO, CHAINF, MINMATCH = args.bedfile, args.from_build, args.to_build, args.chain, args.minmatch

# extract path
PATH = os.path.split(F)[0]
print(PATH, MINMATCH, type(MINMATCH))


# params
LIFTOVER_SRC = "/wynton/group/ahituv/bin/liftOver"

# print message
print("\n\nlifting over", F, "from",  FROM, "to", TO, "in", PATH, "\n\n")


def get_chainFile(from_build, to_build):
    
    # address capitalization issue w goldenpath nomenclature and to-build name
    to_build_caps = to_build[0].upper() + to_build[1:3] + to_build[3].upper() + to_build[4:]
    
    ### get the chainfile between the two builds ###
    
    if to_build == "hs1":  # Hs1 is not yet in wynton group data. Go to my home dir instead. 
        chainPath =f"/wynton/home/ahituv/fongsl/dna/{from_build}"
        
    else:
        chainPath = f"/wynton/group/databases/goldenPath/{from_build}/liftOver" # path to chain file

    
    chainf = os.path.join(chainPath, f"{from_build}To{to_build_caps}.over.chain.gz")
    
    if to_build == "artJam2" and from_build == "hg38":
        # special case for wei's fruit bats
        chainf="/wynton/home/ahituv/fongsl/other_analyses/for-wei_bats/data/chain/hg38.chr1_22.ArtJamW.chain.final"
    
    return chainf
    

def filesToWrite(sid, to_build, path):
    """
    Return file names for liftOver and not liftOver file to be written
    
    input
        sid - (str) sample id
        to_build - (str) build to liftOver to
        
    output
        lifted - (str) file w/ liftOver results
        notlifted - (str) file w/ no liftOver results
    """

    lifted = os.path.join(path, f"{sid}.liftOver.to.{to_build}.bed") # name the liftover file
    notlifted = os.path.join(path, f"{sid}.notlifted.to.{to_build}.bed") # name the notlifted file
    
    return lifted, notlifted
 
    
def sortBed(bedfile, sid, path):
    """
    return a sorted, temporary bed file for lifting over.
    
    input
        bed file (.bed)
        sample id (str)
    method
        1. make tempbed file
        2. compile command line sort command
        3. call sort command
        4. return sorted temp.bed file. 
        
    """
    #1
    tempbed = os.path.join(path, f"temp_{sid}.bed") # sorted temporary bed file

    # [[chr start end enh_id sample_id]] and sort by coordinates
    #2
    cmd = f"sort -k1,1 -k2,2 -k3,3 {bedfile} > {tempbed}"
    
    if os.path.exists(tempbed) is False:
        print("Sorting .bed", tempbed)
    
        #3
        subprocess.call(cmd, shell=True)
    else:
        print("sorted bed already")
    
    #4
    return tempbed
    


def main(argv):
    
    """
    Return liftover and not liftedover .bed files
    
    input 
        .bed file (str) to be lifted over
        path (str) - directory to write liftedover files
        from_build (str) - current build of the input.bed file
        to_build (str) - genome build to liftOver to
    
    method
        1. get sample id
        2. sort .bed file in command line -> write sorted file as temp_file.bed
        3. get chain file for liftOver
        4. get files to write for liftOver, not liftOver results
        5. do the liftover
        6. if liftOver worked, remove the temporary bed file. 
        
    output
        liftedOver file and notliftedOver file 
        
    """
    #1 
    sid = os.path.split(F)[1].strip(".bed")

    #2
    sorted_bed = sortBed(F, sid, PATH)
    
    #3
    if CHAINF == None:
        CHAIN = get_chainFile(FROM, TO)
    else:
        CHAIN=CHAINF

    #4
    lifted, notlifted = filesToWrite(sid, TO, PATH)
    
    #5
    cmd = " ".join([LIFTOVER_SRC, 
                    sorted_bed, 
                    CHAIN,  
                    lifted,  
                    notlifted, 
                    f"-minMatch={MINMATCH}"
                   ])
   
    # if file doesn't exist
    if os.path.exists(lifted) is False:
        print("liftingOver", sid, "\n\n", cmd, "\n\n")
    
        subprocess.call(cmd, shell=True)
    
        print("done lifting")

    # if file exists, but is empty    
    elif os.stat(lifted).st_size == 0:
        print("found empty file\n\nliftingOver", sid, "\n\n", cmd, "\n\n")
    
        subprocess.call(cmd, shell=True)
    
        print("done lifting")
        
    else:
        print("lifted this already?\n\n", lifted)

    #6
    if os.path.exists(lifted) is True:
        if os.path.getsize(lifted) >0:
            os.remove(sorted_bed)
            print("cleaned up temp file")
    

if __name__ == "__main__":
    main(sys.argv[1:])

    
