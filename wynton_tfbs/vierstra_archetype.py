import glob
import os, sys

# input file to intersect against. 
A_BEDFILE = sys.argv[1]
    #e.g. A_BEDFILE = "/wynton/group/ahituv/biomarin/library_1/Design/unq_tiles.lib1.bed"

# vierstra chromosome TFBS archetype annotation files. Large. Zipped. 
PATH = "/wynton/group/ahituv/data/tfbs_motif/vierstra_archetypes_2020/"
vierstra_chrs = glob.glob(os.path.join(PATH, "chr*.bed.gz"))
print(vierstra_chrs)

# make file to write results to 
OUTFILE = os.path.splitext(A_BEDFILE)[0] + ".x.vierstra.archetypes.bed"

# touch outfile
if os.path.exists(OUTFILE) is False:
    os.system(f"touch {OUTFILE}")

# per zipped chromosome archetype file, perform intersection
for B_BEDFILE in vierstra_chrs:
    cmd = " ".join([
        "bedtools intersect -a",
        A_BEDFILE,
        "-b", 
        B_BEDFILE,
        "-wao >>",  # write b info, append to outfile
        OUTFILE
        ])

    print(cmd)

    os.system(cmd)
