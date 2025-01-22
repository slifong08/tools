import argparse
from Bio.SeqUtils import GC as gc_fraction
import glob
import pandas as pd
import os, sys
from statsmodels.stats.multitest import fdrcorrection


###
#   arguments
###
parser = argparse.ArgumentParser(description = "input fasta or bed file, run fimo, get ")

parser.add_argument("input_file", 
                    type=str, 
                    help = ".bed or .fa file for fimo TFBS")

parser.add_argument("-c", "--config", 
                    default=None, type=str, required=False, 
                    help = "config file (optional)")

parser.add_argument("-g", "--genome_build", 
                    default="hg38", type=str, required=False,  
                    help = "genome build (default = hg38)")
# parse arguments
args = parser.parse_args()

# arg variables
INPUT, CONFIG, GENOME_BUILD = args.input_file, args.config, args.genome_build

# str splitting on variables
OUTDIR, SAMPLE_NAME = os.path.split(INPUT)
SAMPLE_ID= os.path.splitext(SAMPLE_NAME)[0]
FIMO_RESULT_DIR = os.path.join(OUTDIR, f"fimo")
print("OUTPUT_DIR", FIMO_RESULT_DIR)

# FILE CONSTANTS
FIMO_SRC = "/wynton/group/ahituv/bin/meme-5.5.5/src/fimo"
MEME = "/wynton/group/ahituv/data/tfbs_motif/jaspar/JASPAR2022_CORE_non-redundant_pfms_meme.txt"


## MK FIMO RESULTS DIR
#if os.path.exists(os.path.join(OUTDIR, f"fimo")) is False and os.path.exists(FIMO_RESULT_DIR) is False:
#    os.mkdir(os.path.join(OUTDIR, f"fimo"))

### 
# FUNCTIONS
###

def bedtoolsGetFasta(bed_file, genome_build):
    """ 
    convert .bed to .fa for specific genome_build
    
    input args
        bed_file (str)
            path to bed file
            
        genome_build (str) 
            genome_build for reference fa
        
    requirements 
        bedtools in local bash environment
        
    method
        1. make handle for fasta file
        2. get ref genome fasta file
        3. run bedtools getfasta command
            if .fa file does not exist and reference fa file is available. 
            else throws error
    
    return
        fa_file (str) - fasta file from bed file. 
    """
    
    # 1
    fa_file = os.path.splitext(bed_file)[0] + ".fa"

    # 2 get genome reference .fa file
    if genome_build in ["hg38", "hg19", "mm10", "mm9"]:

        #ref_fasta = f"/wynton/group/databases/goldenPath/{genome_build}/bigZips/{genome_build}.fa.gz"

        # ref fasta must be unzipped
        ref_fasta = f"/wynton/group/ahituv/ref/{genome_build}/Sequence/WholeGenomeFasta/genome.fa"
    else:
        print("no reference genome build fasta file! update bedtoolsGetFasta function")
        ref_fasta = None
    
    # 3
    if os.path.exists(fa_file) is False and ref_fasta is not None:
        
        # use bedtools getFasta to go from .bed -> .fa
        cmd = " ".join([
            "bedtools getfasta",
            "-fi",
            ref_fasta,
            "-bed",
            bed_file,
            ">",
            fa_file
        ])
        print("running bedtools getfasta for", bed_file )
        os.system(cmd)
        
    elif ref_fasta is not None:
        print("bed is already fasta")
        
    return fa_file


def runFimo(fa, sample_id, fimo_src, meme):
    
    """
    command to run fimo
    
    input args
        fa (str) 
            fasta file w full path
        sample_id (str)
            name of the sample (filename)
        fimo_src (str)
            fimo src code to run fimo search
        meme (str)
            meme files for motifs
            
    method
        1. get fa path, make output_dir
        2. chdir to dir where .fa is
        3. compile fimo command 
               --no-qvalue: no BH q-value pred. Do this later
               -- text: write only .tsv
        4. run command through the commandline
            only if not run before
    
    return
        output_dir (str)
            name of directory where fimo outputs live. 
    """
    
    #1
    dir_, fa_name = os.path.split(fa) 
    output_dir = os.path.join(dir_, f"fimo")
    output_file= os.path.join(output_dir, "fimo", "fimo.txt")
    
    #2
    os.chdir(dir_)

    #3
    cmd = " ".join([
                    fimo_src,
                    "--no-qvalue --skip-matched-sequence --verbosity 1",  # do not compute FDR, output tsv, do not print matches
                    meme,
                    fa, 
                    ">",
                    output_file,

                ])
    
    #4

    print(cmd)
    os.system(cmd)
    

    return output_dir


def fimo2bed(fimo_outdir, sample_id):
    

    
    """
    format fimo outputs after run
    
    inputs
        fimo_outdir (str) 
            directory w fimo data
        sample_id (str) 
            name of file
    
    method
        1. make an outfile, check if file exists already.
        2. if not, process raw fimo output file. 
        3. open raw file as df
        4. TF-name col from str split
        5. element genome coordinate cols
        6. motif relative coor -> motif genome coor cols
        7. GC, len cols 
        8. rearrange col order
        9. 5% FDR correction for p-values per motif. 
        10. write df
        11. throw message if file already exists. 
    
    return
        None
    """
    
    #1
    outfile = os.path.join(fimo_outdir, f"{sample_id}.post.run.motifs.bed")
    
    # 2
    if os.path.exists(outfile) is False:
        # 3
        df = pd.read_csv(f"{fimo_outdir}/fimo.txt", sep = '\t')

        # 4 TF name
        df["tf"] = df["motif_alt_id"].apply(lambda x: x.split("_")[0])

        # 5 genome coor
        df["#chr_element"] = df["sequence_name"].apply(lambda x: x.split(":")[0])
        df["start_element"] = df["sequence_name"].apply(lambda x: (x.split(":")[1]).split("-")[0])
        df["end_element"] = df["sequence_name"].apply(lambda x: (x.split(":")[1]).split("-")[1])
        
        # 6 motif relative coor -> motif genome coor 
        df["start_motif"] = df["start_element"].astype(int) + df["start"].astype(int)
        df["end_motif"] = df["end_element"].astype(int) + df["stop"].astype(int)

        # 7 GC, len features of motif
        df["motif_GC"] = df["matched sequence"].apply(lambda x: gc_fraction(x))
        df["motif_len"] = df["matched sequence"].apply(lambda x: len(x))

        # 8 rearrange columns
        df = df[['#chr_element',"start_motif", "end_motif", 
                 'start_element', 'end_element', 
                 'tf', 'motif_alt_id',
                 'sequence_name', 'start', 'stop', 'strand', 'score',
                 'p-value','motif_GC', 'motif_len', 'matched sequence'
                ]].drop_duplicates()
        
        #9 perform per motif FDR correction
        memes = set(df['motif_alt_id'])
        
        fdr_dict={}
        for meme in memes:
            test = df.loc[df["motif_alt_id"] == meme]
            test["FDR_bool"], test["FDR_p"] = fdrcorrection(test["p-value"])
            fdr_dict[meme] = test
        
        fdr = pd.concat(fdr_dict.values())
       
        # 10 write
        fdr.to_csv(outfile, sep = "\t", header = False, index = False)
        print("made!", outfile)
    # 11
    else:
        print("already made?", outfile)

    return outfile
    
def writeConfig(config_name, fimo_src, meme, results_dir):
    
    """ 
    import config_readwrite function, 
    write fimo, meme, output dir to config file
    """
    
    sys.path.append("/wynton/group/ahituv/fongsl/tools/py_")
    import config_readwrite as crw
    
    config, cfn = crw.read(config_name) 

    # write to config
    section = f"FIMO"
    crw.check(config, section)

    # fimo src code path
    config[section]["src"] = fimo_src

    # JASPAR meme file path 
    config[section]["meme"] = meme
    

    # directory to FIMO results
    config[section]["path"] = results_dir

    # save
    crw.write(config, cfn)
    print("wrote config")

    
###
#   main
###

def main(argv):
    """ 
    run fimo, write to config, process fimo outputs
    
    method 
        1. if config is not None, write fimo variables to config
        2. if .bed file input, convert to fasta w bedtools getfasta
        3. run fimo
        4. process fimo outputs
    """
    
    # 1
    if CONFIG is not None:
        writeConfig(CONFIG, FIMO_SRC, MEME, FIMO_RESULT_DIR)
    
    # 2
    if ".bed" in INPUT:
        FA = bedtoolsGetFasta(INPUT, GENOME_BUILD)
    elif ".fa" in INPUT:
        FA = INPUT
    else:
        print("input is not .bed or .fa")
    # 3   
    FIMO_OUTDIR = runFimo(FA, SAMPLE_ID, FIMO_SRC, MEME) # run fimo
    
    # 4
    #xout = fimo2bed(FIMO_OUTDIR, SAMPLE_ID) # fimo text to bed

if __name__ == "__main__":
    main(sys.argv[1:])
