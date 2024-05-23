#!/bin/bash

#SBATCH --job-name=test_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --time=0-01:00:00
#SBATCH --mem=256GB
#SBATCH --error=/data/biol-micro-genomics/kell7366/plink/test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/plink/test_output.txt

#import modules
module purge
module load Anaconda3
source activate $DATA/myenv
#module load BWA/0.7.17-GCCcore-11.2.0
#module load PLINK/2.00a2.3_x86_64
#module load SAMtools/1.16.1-GCC-11.3.0
#module load BCFtools/1.14-GCC-11.2.0
module load R/4.2.2-foss-2022a

rgenome="GCA_000465235.1_ASM46523v1_genomic.fna"
bam_output="/data/biol-micro-genomics/kell7366/plink/coli_bam"

#mkdir -p "$bam_output"
#
## Index the reference genome
#bwa index "$rgenome"
#samtools faidx "$rgenome"

## Align each isolate to the reference genome, add unique read group, and convert to BAM
#for isolate in jejuni/*.fas; do
#  bam_name=$(basename "$isolate" .fas)
#  read_group="@RG\tID:$bam_name\tSM:$bam_name\tPL:ILLUMINA" # Assuming Illumina platform
#
#  if [ -e "$bam_output/${bam_name}.bam" ]; then
#    echo "File ${bam_name}.bam already exists. Skipping."
#  else
#    echo "Processing $bam_name"
#    bwa mem -R "$read_group" "$rgenome" "$isolate" | \
#    samtools view -Sb - > "$bam_output/${bam_name}.bam"
#  fi
#done
#
## Merge BAM files while retaining read group information
#samtools merge -r "$bam_output/merged.bam" "$bam_output"/*.bam
#
## Sort and index the merged BAM file
#samtools sort "$bam_output/merged.bam" -o "$bam_output/merged.sorted.bam"
#samtools index "$bam_output/merged.sorted.bam"
#
## Variant calling using samtools and bcftools on the merged and sorted BAM file
#bcftools mpileup -Ou -f "$rgenome" "$bam_output/merged.sorted.bam" | \
#bcftools call -mv -Ob -o "$bam_output/variants.bcf"
#
#bcftools view "$bam_output/variants.bcf" -o "$bam_output/variants.vcf"

## Assuming variants.vcf is already generated
#
## Count the number of SNPs (non-header lines)
#total_snps=$(grep -vc "^#" "$bam_output/variants.vcf")
#
## Calculate 1/10th of the total SNPs
#target_snps=$((total_snps / 3))
#
## Randomly select approximately 1/10th of the SNPs
## Keep the VCF header
#grep "^#" "$bam_output/variants.vcf" > "$bam_output/variants_reduced.vcf"
#
## Randomly select SNPs and add them to the reduced VCF file
#grep -v "^#" "$bam_output/variants.vcf" | shuf | head -n "$target_snps" >> "$bam_output/variants_reduced.vcf"
#
## PLINK command to calculate LD (assuming VCF file is generated)
#./plink_linux_x86_64_20231211/plink --double-id --allow-extra-chr --vcf "$bam_output/variants_reduced.vcf"  --make-bed --out temp
#
#./plink_linux_x86_64_20231211/plink --double-id --allow-extra-chr --bfile temp --ld-window 2000000 --threads 12 --r2 --ld-window-kb 1600 --ld-window-r2 0 --out ld_output

#awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' ld_output.ld > ld_output_tab_delimited.ld
#
#rm ld_output.ld

#awk 'BEGIN {srand()} NR==1 || (!/^($|[:space:]*#)/ && rand() <= 0.01) { print $0 }' ld_output_tab_delimited.ld > ld_output_sampled.ld
#
## Define the directory containing ld_chunk files
#ld_chunk_dir="split"
#
#mkdir -p "$ld_chunk_dir"
#
#split -l 10000000 ld_output_tab_delimited.ld ./split/ld_chunk_
#
#rm ld_output_tab_delimited.ld

### Define an array of specific files to process
#declare -a files_to_process=("split/ld_chunk_bn") # Add your file names here
#
## Process each specified file
#for file in "${files_to_process[@]}"; do
##for file in "${ld_chunk_dir}/"/*; do
#    echo "Processing file: $file"
#    awk -F'\t' -v OFS='\t' '
#    NR==FNR {
#        # Reading gene_target.tsv and storing data in an array
#        gene_targets[$1] = $3 OFS $4 OFS $6 OFS $7
#        next
#    }
#    {
#        # Processing each line of ld_chunk files
#        for (link_name in gene_targets) {
#            split(gene_targets[link_name], genes, OFS)
#            if (($2 >= genes[1] && $2 <= genes[2]) && ($5 >= genes[3] && $5 <= genes[4])) {
#                print $1, $2, ".", $3, $5, ".", $7 >> "link_"link_name".tsv"
#            }
#        }
#    }' gene_target.tsv "$file"
#done
#
#echo "Processing complete."

Rscript plotting_lddecay.r