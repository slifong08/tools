#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l h_vmem=16G      # memory 
#$ -l mem_free=16G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=16G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=12:15:00   # job requires up to 24 hours of runtime
##$ -t 1-10           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to
#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".output.txt
##$ -t 1-2:1        # job range:step size
##$ -tc 1             # n jobs to run at once

## If you array jobs (option -t), this script will run T times, once per task.
## For each run, $SGE_TASK_ID is set to the corresponding task index (here 1-10).
## To configure different parameters for each task index, one can use a Bash 
## array to map from the task index to a parameter string.

## All possible parameters
# params=(1bac 2xyz 3ijk 4abc 5def 6ghi 7jkl 8mno 9pqr 10stu)

## Select the parameter for the current task index
## Arrays are indexed from 0, so we subtract one from the task index
# param="${params[$((SGE_TASK_ID - 1))]}"

date
hostname

#ml load cuda/11.5
micromamba activate /wynton/home/ahituv/fongsl/micromamba/envs/mamba


PORT=7778
echo To open a tunnel from local machine, 
echo Execute in a new terminal window:
echo "ssh -L 9999:$HOSTNAME:$PORT $(whoami)@log1.wynton.ucsf.edu"
printf '=%.0s' {1..80}
echo


jupyter notebook --no-browser --ip='*' --port=${PORT}

# $1 is .bed file to intersect

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
