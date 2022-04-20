#!/bin/bash
#SBATCH -J "nf-scarches"
#SBATCH --output=/data/sennis/AML/logs/nf-scarches_test.out
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 8

nextflow='/home/sennis/nextflow'
cd /data/sennis/AML/scarches_nf
module load singularity
$nextflow run main.nf -profile singularity -with-trace trace.txt -with-dag flowchart.png
