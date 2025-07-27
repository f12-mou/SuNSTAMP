#!/bin/bash
#SBATCH --job-name=faria_ssp-conda
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=23:30:00
#SBATCH --chdir=/mount/samee/faria-suns/SunSTAMP/replicate_8_1/src
#SBATCH --error %x.%A_%a.joberr
#SBATCH --output %x.%A_%a.jobout
#SBATCH --array=0-20

### enable 'module' commands
source /etc/profile.d/lmod.sh

### Allow common binaries from conda
module load anaconda3/2020.11
source /home/samee/.bashrc

pwd
conda activate m_faria
/mount/samee/.conda/envs/m_faria/bin/python3 nmf_final_2.py
pwd
