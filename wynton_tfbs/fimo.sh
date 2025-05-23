#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=35G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=35G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=12:15:00   # job requires up to 24 hours of runtime
##$ -t 1-10           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to
#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".output.txt


# ONLY FOR HG38

echo micromamba activate $HOME/micromamba/envs/mm2025
micromamba activate $HOME/micromamba/envs/mm2025

echo python3 $HOME/tools/wynton_tfbs/fimo.py "$1" 
python3 $HOME/tools/wynton_tfbs/fimo.py "$1" 
# $SGE_TASK_ID = chrnum

# $1 = hg38 fasta file!


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"