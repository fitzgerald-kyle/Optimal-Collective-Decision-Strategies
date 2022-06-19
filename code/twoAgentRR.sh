#!/bin/bash

#SBATCH --qos=blanca-appm
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=twoAgentRR
#SBATCH --output=twoAgentRR_%j.out
#SBATCH --mail-user=kyfi8465@colorado.edu
#SBATCH --mail-type=ALL

module purge
module load matlab
matlab -nodesktop -nondisplay -r 'clear; twoAgentRR;'
