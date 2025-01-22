#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=2G     # job requires up to 1 GiB of RAM per slot
#$ -l h_rt=08:00:00   # job requires up to 24 hours of runtime
#$ -r y               # if job crashes, it should be restarted


cd /Wynton/group/ahituv/
chmod 777 -R ./ALU
chmod 777 -R ./nadja
chmod 777 -R ./liver_nullomer_JJ_analysis
chmod 777 -R ./MPRAbase
chmod 777 -R ./MPRAflow_old
chmod 777 -R ./MPRAhub-LZ
chmod 777 -R ./MultiSeq
chmod 777 -R ./Ofer_ATAC
chmod 777 -R ./Ofer_RNA
chmod 777 -R ./Rib
chmod 777 -R ./SCN2A_RNASEQ