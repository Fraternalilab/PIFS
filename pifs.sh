#!/bin/bash
# 
#SBATCH -J popscomp                                                             # Job name
#SBATCH -p long                                                                 # partition (queue)
#SBATCH -N 1                                                                    # number of nodes 
#SBATCH -n 1                                                                    # number of cores 
#SBATCH -o slurm.%j.out                                                         # STDOUT 
#SBATCH -e slurm.%j.err                                                         # STDERR 

#############################################################################
# Run PIFS over a list of PDBs. Large structures (>20 chains) were ignored. #
#############################################################################

python -u pifs.py --structList $file --inDir ../pdb/20180101/ --networkOutDir $networkDir --popscompOutDir ../pdbInterfaceScreen/popscomp/ --pifsDir $pifsDir --pifsTag $pifsTag --popsDir ../ --popscompDir ../POPSCOMP/ --errorFile /data/home/josephn/pifs/`echo $pifsTag`.err --distCutoff 6.0
