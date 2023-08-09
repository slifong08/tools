import argparse
import os, sys
import subprocess

###
# arguments
###

arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("-b","--bedfile",  type=str,help='bed file w/ full path')
arg_parser.add_argument("-f","--from_build", type=str, help='genome build of the original file')
arg_parser.add_argument("-t","--to_build", type=str, help='genome build to liftover to')

args = arg_parser.parse_args()

F = args.bedfile
from_build= args.from_build
to_build= args.to_build

path = "/".join(F.split("/")[:-1])

print("lifting over", F, "from, to", from_build, to_build, "in", path)

def get_chainFile(from_build, to_build):
    
     ### get the chainfile between the two builds ###

    chainPath = "/dors/capra_lab/data/ucsc/liftOver/" # path to chain file
    chainf = os.path.join(chainPath, f"{from_build}To{to_build}.over.chain.gz")
    
    return chainf
    

def filesToWrite(sid, to_build):
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
    
def sortBed(bedfile, sid):
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

    print("Sorting .bed", tempbed)
    
    #3
    subprocess.call(cmd, shell=True)
    
    #4
    return tempbed
    
def liftover(bedfile, path, from_build, to_build): # bedfile with full path
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
    sid = (bedfile.split("/")[-1]).split(".")[0] # get the sample ID

    #2
    sorted_bed = sortBed(bedfile, sid)
    
    #3
    chainf = get_chainFile(from_build, to_build)
    
    #4
    lifted, notlifted = filesToWrite(sid, to_build)
    
    #5
    cmd = f"liftOver {sorted_bed} {chainf} {lifted} {notlifted}"
    print("liftingOver", sid, "\n\n", cmd, "\n\n")
    
    subprocess.call(cmd, shell=True)
    
    print("done lifting")

    #6
    if os.path.getsize(lifted) >0:
        os.remove(sorted_bed)
        print("cleaned up temp file")

    return lifted

def main(argv):
    
    liftover(F, path, from_build, to_build)
#%%

if __name__ == "__main__":
    main(sys.argv[1:])
