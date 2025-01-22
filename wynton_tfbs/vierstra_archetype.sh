#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l h_vmem=50G      # memory 
#$ -l mem_free=50G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=50G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=12:15:00   # job requires up to 24 hours of runtime
##$ -t 1-10           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to
#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".output.txt
##$ -t 1-2:1        # job range:step size
##$ -tc 1             # n jobs to run at once


date
hostname
ml load CBI  bedtools2/2.26.0
micromamba activate /wynton/home/ahituv/fongsl/.conda/envs/mamba # environment

echo "$HOME"/vierstra_archetype.sh

python3 "$HOME"/vierstra_archetype.py $1

# $1 is .bed file to intersect

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
