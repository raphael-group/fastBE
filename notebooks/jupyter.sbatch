#!/bin/bash
#SBATCH --mincpus 32
#SBATCH --mem 216G
#SBATCH --time 6-23:59:00
#SBATCH --job-name jupyter
#SBATCH --output jupyter_logs/jupyter-notebook-%J.log
# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i10001-20000 -n1)
node=$(hostname -s)
user=$(whoami)
cluster="soak"

# print tunneling instructions jupyter-log
echo -e "
MacOS or linux terminal command to create your ssh tunnel
ssh -N -L ${port}:${node}:${port} ${user}@${cluster}.cs.princeton.edu
Windows MobaXterm info
Forwarded port:same as remote port
Remote server: ${node}
Remote port: ${port}
SSH server: ${cluster}.cs.princeton.edu
SSH login: $user
SSH port: 22
Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"

# load modules or conda environments here
module load gurobi
source /n/fs/ragr-data/users/schmidt/miniconda3/bin/activate base
conda activate breaked
cd /n/fs/ragr-research/projects/fastBE/

# use one of the below two lines depending on what you want to run
# jupyter lab --no-browser --port=${port} --ip=${node}
jupyter-notebook --no-browser --port=${port} --ip=${node}
