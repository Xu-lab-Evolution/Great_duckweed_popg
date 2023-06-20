#!/bin/bash

#SBATCH -J chr10                  # Job name
#SBATCH -o chr10.%j.out           # Specify stdout output file (%j expands to jobId)
#SBATCH -p smp                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 1                     # Total number of tasks, here explicitly 1
#SBATCH --mem 2G                 # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 10:00:00              # Run time (hh:mm:ss)
 
#SBATCH -A  m2_jgu-evoltroph     # Specify allocation to charge against

#SBATCH --mail-type=END
#SBATCH --mail-user=pduchnbo@uni-mainz.de
 
# Load all necessary modules if needed (these are examples)
# Loading modules in the script ensures a consistent environment.
 
# Launch the executable
exec=/home/pduchnbo/Programs/3P-CLR-master/src/threepclr
srun $exec "input_Chr10.txt" "output_Chr10.txt" 20 100 0.0025 0.0619094587263976,0.048751715849206,0.294460819600533 NA

