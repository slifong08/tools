#!/bin/bash
#SBATCH --mail-user=sarah.l.fong@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1    # comments allowed
#SBATCH --ntasks=1
#SBATCH --time=00:60:00
#SBATCH --output={$6}
#SBATCH --job-name=CON_ACC
#SBATCH --account=capra_lab
#SBATCH --output=%J.stderr


ml load Anaconda3 
source /home/fongsl/.conda/envs/sfenv/bin/activate

echo $SLURM_ARRAY_TASK_ID
echo chr ${1} br ${2} msa ${3} mod ${4} outdir ${5} 
#echo python /dors/capra_lab/users/fongsl/tools/evo/conacc_array.py -i $SLURM_ARRAY_TASK_ID -chr ${1} -br ${2} -g ${3} -msa ${4} -mod ${5} -o ${6}
echo python /dors/capra_lab/users/fongsl/tools/evo/conacc_array.py -i $SLURM_ARRAY_TASK_ID -chr ${1} -br ${2} -g ${3} -msa ${4} -mod ${5} -o ${6}

python /dors/capra_lab/users/fongsl/tools/evo/conacc_array.py -chr ${1} -br ${2} -g ${3} -msa ${4} -mod ${5} -o ${6}