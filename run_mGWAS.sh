#!/bin/bash

#SBATCH --job-name=testrun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=short
#SBATCH --time=0-02:00:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/mGWAS/test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/mGWAS/test_output.txt

bash run_maskfasta.sh -q cassette_genes.fasta -i ../minimap/ref_fas/ref_aureus/positive/ -d 3000 -o . -t 12