import os

def makeWindows(bedfilename, windowsize=270, stepsize=1, dryrun=True):
    
    """return windowed bedfile from input bedfile, stepsize, makewindows software from BedTools.

    ###
    # Requirements
    ###
    1. a bedfile
    2. bedtools software downloaded in linux environment. 

    ###
    # Input
    ###
    - .bed - (str) file of coordinates to window over. Can have a fourth id column
    - windowsize - (int) desired size of window. Default is 270
    - stepsize - (int) desired step/shift size. Default is one bp
    - dryrun - (bool)

    ###
    # Methods
    ###
    1. store outputfile, derived from bedfilename as variable
    2. compile bedtools command
    3. if dryrun=False, run command in commandline.
    4. return windowed.bed file
    
    
    """

    # step 1
    outfile = str(bedfilename).replace(".bed", f".window-{windowsize}.step{stepsize}.bed")

    # step2
    cmd = " ".join([
                    "bedtools makewindows -b", 
                    str(bedfilename), 
                    "-w", 
                    str(windowsize), 
                    "-i srcwinnum" , 
                    "-s", 
                    str(stepsize), 
                    ">", 
                    outfile
                    ])
    print(cmd)

    # step 3
    if dryrun is False:
        os.system(cmd)

    # step 4
    
    return outfile

## params 
OLIGO_LEN=270
STEPSIZE=50
TEST_BED='./test.bed'

# Window the bed file. 
WINDOWBED = makeWindows(bedfilename=TEST_BED, windowsize=OLIGO_LEN, stepsize=STEPSIZE, dryrun=False)

