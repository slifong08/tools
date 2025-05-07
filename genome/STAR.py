# download STARR aligner
import os

# source dir
os.chdir("/wynton/group/ahituv/bin")

STARSRC = "/wynton/group/ahituv/bin/STAR/source/"
GENOME_DIR='/wynton/group/ahituv/data/genome_reference/STAR/hg38/'
FA = "/wynton/group/ahituv/data/dna/hg38/hg38.fa"
GTF = "/wynton/group/ahituv/data/dna/hg38/hg38.knownGene.gtf"


# step 1 - check if STAR has already been downloaded. 
# see https://github.com/alexdobin/STAR/tree/master
if os.path.exists(STARSRC) is False:
    
    # clone git 
    git_cmd= "git clone https://github.com/alexdobin/STAR.git"
    
    os.system(git_cmd)
    
    os.chdir(STARSRC)
    
    os.system("make STAR")


# step 2 - make index

index_genome_cmd = [
    "./STAR --runThreadN 10 --runMode genomeGenerate --genomeDir", 
    GENOME_DIR, 
    "--genomeFastaFiles", 
    FA, 
    "--sjdbGTFfile", 
    GTF
]

os.system(" ".join(index_genome_cmd))
