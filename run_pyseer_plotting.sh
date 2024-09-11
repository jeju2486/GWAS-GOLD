#!/bin/bash

#module load
module purge
module load Anaconda3/2024.02-1
source activate $DATA/python3_8_18

pyseer_github_scripts="/data/biol-micro-genomics/kell7366/pyseer/pyseer_codes/scripts"

reference_seq="/data/biol-micro-genomics/kell7366/sccmec_ncbi_dataset/filtered_data/ref_fas/ref_aureus/all/GCF_000144955.fas"

python "$pyseer_github_scripts"/../phandango_mapper-runner.py "$result_dir"/sccmec_kmers.filtered.txt "$reference_seq" "$result_dir"/saureus_kmers_minimum_GCF_000144955.plot
python "$pyseer_github_scripts"/../phandango_mapper-runner.py "$result_dir"/high_outliers_kmers.txt "$reference_seq" "$result_dir"/saureus_kmers_high_GCF_000144955.plot

# This will iterate down the list of annotations, annotating the k-mers which haven뭪 already been mapped to a previous annotation
python "$pyseer_github_scripts"/../annotate_hits_pyseer-runner.py "$result_dir"/filtered_kmers_normal_and_low.txt reference_aureus.txt "$result_dir"/low_annotated_kmers.txt
python "$pyseer_github_scripts"/../annotate_hits_pyseer-runner.py "$result_dir"/high_outliers_kmers.txt reference_aureus.txt "$result_dir"/high_annotated_kmers.txt

# summarise annotations to create a plot of significant genes, only use genes k-mers are actually in
python "$pyseer_github_scripts"/summarise_annotations.py "$result_dir"/low_annotated_kmers.txt > "$result_dir"/low_gene_hits.txt
python "$pyseer_github_scripts"/summarise_annotations.py "$result_dir"/high_annotated_kmers.txt > "$result_dir"/high_gene_hits.txt

module purge
module load R/4.2.2-foss-2022a

Rscript plotting_annotation.R "$result_dir"/low_gene_hits.txt low_annotation
Rscript plotting_annotation.R "$result_dir"/high_gene_hits.txt high_annotation