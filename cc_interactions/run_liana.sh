#!/bin/bash
#SBATCH -J "liana"
#SBATCH --output=/data/sennis/AML/logs/liana.out
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8

cd /data/sennis/AML/bone_marrow_analysis
module load singularity
singularity exec -B /data cc_interactions/cci_image.simg cc_interactions/run_liana.R /data/sennis/AML/bone_marrow_analysis/data/healthy_sub.Rds
singularity exec -B /data cc_interactions/cci_image.simg cc_interactions/run_liana.R /data/sennis/AML/bone_marrow_analysis/data/d0_sub.Rds
