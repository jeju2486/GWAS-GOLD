#!/bin/bash

#SBATCH --job-name=ctmc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=short
#SBATCH --time=0-02:00:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/mGWAS_staphy/mGWAS_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/mGWAS_staphy/mGWAS_output.txt

# Load necessary modules
module purge
module load minimap2/2.24-GCCcore-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0
module load BWA/0.7.17-GCCcore-11.2.0
module load Anaconda3/2023.09-0
source activate $DATA/python3_12_2

#directory calling
query_seq="/data/biol-micro-genomics/kell7366/maskGWAS/clinical_blaCTX-Ms.fasta"
input_dir="/data/biol-micro-genomics/kell7366/maskGWAS/input_fas"
maskfasta_output_dir="/data/biol-micro-genomics/kell7366/maskGWAS/output"

echo "CPU numbers: $SLURM_CPUS_PER_TASK"

#Run the code
#bash run_maskfasta.sh -q "$query_seq" -i "$input_dir" -d 5151 -o "$maskfasta_output_dir" -t $SLURM_CPUS_PER_TASK

echo "Running the pyseer"

#Reload necessary modules
conda deactivate
source activate $DATA/python3_8_18

# Define paths for Pyseer
treefile="/data/biol-micro-genomics/kell7366/maskGWAS/phylogroup_F_core_tree.treefile"
phenotype="/data/biol-micro-genomics/kell7366/maskGWAS/phylogroup_F_phenotype.tsv" #It should be tab delimited and two column file
pyseer_output_dir="/data/biol-micro-genomics/kell7366/maskGWAS/pyseer_output"
pyseer_script="/data/biol-micro-genomics/kell7366/pyseer/pyseer_codes/scripts"

## Run the Pyseer script
#bash run_pyseer.sh \
#   -t "$treefile" \
#   -p "$phenotype" \
#   -i "$maskfasta_output_dir/unitig_output" \
#   -o "$pyseer_output_dir" \
#   -P "ecoli" \
#   -T "$SLURM_CPUS_PER_TASK" \
#   -s "$pyseer_script"

echo "Generating the phandango plots"

reference_seq="/data/biol-micro-genomics/kell7366/maskGWAS/input_fas/20001_562_20504.fasta"

python "$pyseer_script"/../phandango_mapper-runner.py "$pyseer_output_dir"/ecoli_kmers.filtered.txt "$reference_seq" "$pyseer_output_dir"/ecoli_kmers_minimum_F.plot
python "$pyseer_script"/../phandango_mapper-runner.py "$pyseer_output_dir"/high_outliers_kmers.txt "$reference_seq" "$pyseer_output_dir"/ecoli_kmers_high_F.plot

## This will iterate down the list of annotations, annotating the k-mers which haven? already been mapped to a previous annotation
python "$pyseer_script"/../annotate_hits_pyseer-runner.py "$pyseer_output_dir"/filtered_kmers_normal_and_low.txt "$pyseer_output_dir"/reference.tsv "$pyseer_output_dir"/low_annotated_kmers.txt
python "$pyseer_script"/../annotate_hits_pyseer-runner.py "$pyseer_output_dir"/high_outliers_kmers.txt "$pyseer_output_dir"/reference.tsv "$pyseer_output_dir"/high_annotated_kmers.txt

# summarise annotations to create a plot of significant genes, only use genes k-mers are actually in
python "$pyseer_script"/summarise_annotations.py "$pyseer_output_dir"/low_annotated_kmers.txt > "$pyseer_output_dir"/low_gene_hits.txt
python "$pyseer_script"/summarise_annotations.py "$pyseer_output_dir"/high_annotated_kmers.txt > "$pyseer_output_dir"/high_gene_hits.txt

module purge
module load R/4.4.0-gfbf-2023a

Rscript plotting_annotation.R "$pyseer_output_dir"/low_gene_hits.txt low_annotation
Rscript plotting_annotation.R "$pyseer_output_dir"/high_gene_hits.txt high_annotation

echo "All done!"