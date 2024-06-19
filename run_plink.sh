#!/bin/bash

# Function to display usage message
usage() {
  echo "Usage: $0 -i <isolate_dir> -o <bam_output>"
  exit 1
}

# Parse input arguments
while getopts "i:o:" opt; do
  case ${opt} in
    i )
      isolate_dir=$OPTARG
      ;;
    o )
      bam_output=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

# Check if the required arguments are provided
if [ -z "$isolate_dir" ] || [ -z "$bam_output" ]; then
  usage
fi

# Load necessary modules
module purge
module load BWA/0.7.17-GCCcore-11.2.0
module load SAMtools/1.16.1-GCC-11.3.0
module load BCFtools/1.14-GCC-11.2.0

# Choose the first .fna, .fasta, or .fas file in the isolate directory as the reference genome
rgenome=$(ls "$isolate_dir"/*.{fna,fasta,fas} 2>/dev/null | head -n 1)

# Check if a reference genome was found
if [ -z "$rgenome" ]; then
  echo "No .fna, .fasta, or .fas files found in the isolate directory."
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$bam_output"

# Index the reference genome
bwa index "$rgenome"
samtools faidx "$rgenome"

# Align each isolate to the reference genome, add unique read group, and convert to BAM
for isolate in "$isolate_dir"/*.{fna,fasta,fas}; do
  # Check if the file exists to avoid errors when there are no matching files
  if [ -e "$isolate" ]; then
    bam_name=$(basename "$isolate")
    bam_name="${bam_name%.*}"  # Remove file extension
    read_group="@RG\tID:$bam_name\tSM:$bam_name\tPL:ILLUMINA" # Assuming Illumina platform

    if [ -e "$bam_output/${bam_name}.bam" ]; then
      echo "File ${bam_name}.bam already exists. Skipping."
    else
      echo "Processing $bam_name"
      bwa mem -R "$read_group" "$rgenome" "$isolate" | \
      samtools view -Sb - > "$bam_output/${bam_name}.bam"
    fi
  fi
done

# Merge BAM files while retaining read group information
samtools merge -r "$bam_output/merged.bam" "$bam_output"/*.bam

# Sort and index the merged BAM file
samtools sort "$bam_output/merged.bam" -o "$bam_output/merged.sorted.bam"
samtools index "$bam_output/merged.sorted.bam"

# Variant calling using samtools and bcftools on the merged and sorted BAM file
bcftools mpileup -Ou -f "$rgenome" "$bam_output/merged.sorted.bam" | \
bcftools call -mv -Ob -o "$bam_output/variants.bcf"

bcftools view "$bam_output/variants.bcf" -o "$bam_output/variants.vcf"

# Assuming variants.vcf is already generated

# Count the number of SNPs (non-header lines)
total_snps=$(grep -vc "^#" "$bam_output/variants.vcf")

# Calculate 1/10th of the total SNPs
target_snps=$((total_snps / 10))

# Randomly select approximately 1/10th of the SNPs
# Keep the VCF header
grep "^#" "$bam_output/variants.vcf" > "$bam_output/variants_reduced.vcf"

# Randomly select SNPs and add them to the reduced VCF file size
grep -v "^#" "$bam_output/variants.vcf" | shuf | head -n "$target_snps" >> "$bam_output/variants_reduced.vcf"

# PLINK command to calculate LD (assuming VCF file is generated)
if [ -e "$bam_output"/ld_output_tab_delimited.ld ]; then  
    
    echo "ld_output already exists. Skipping PLINK commands."
    
  else

  ./plink_linux_x86_64_20231211/plink --double-id --allow-extra-chr --vcf "$bam_output/variants_reduced.vcf" --make-bed --out temp
  ./plink_linux_x86_64_20231211/plink --double-id --allow-extra-chr --bfile temp --ld-window 200000 --threads 1 --r2 --ld-window-kb 2900 --ld-window-r2 0 --out "$bam_output"/ld_output
  
  # Make file tab-delimited for better processing after
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' "$bam_output"/ld_output.ld > "$bam_output"/ld_output_tab_delimited.ld
  
  rm temp*
  rm "$bam_output"/ld_output.ld

fi

# For reducing the time of plotting let sampling 10% of total sample
echo "Sampling for LD graph starts..."

awk 'BEGIN {srand()} NR==1 || (!/^($|[:space:]*#)/ && rand() <= 0.05) { print $0 }' "$bam_output"/ld_output_tab_delimited.ld > "$bam_output"/ld_output_sampled.ld

for index_file in "$isolate_dir"/*.{amb,ann,bwt,fai,pac,sa}; do
  mv $index_file "$bam_output"
done

echo "Processing complete."

echo "Plotting start..."

# Load R
module purge
module load R/4.2.2-foss-2022a

Rscript plotting_lddecay.r -i "$bam_output"/ld_output_sampled.ld -o "$bam_output"