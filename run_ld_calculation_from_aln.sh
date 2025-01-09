#!/usr/bin/env bash
set -euo pipefail

################################################################################
# run_ld_calculation.sh
#
# 1) Converts Pirate's core alignment to a multi-sample VCF using snp-sites
#    (if not already done).
# 2) Randomly selects 5% of the variants from that VCF to reduce computational
#    burden on PLINK.
# 3) Converts the subsampled VCF to PLINK format.
# 4) Performs pairwise LD analysis in PLINK.
# 5) Outputs LD results (and a further 5% subsample of the LD pairs for plotting,
#    if you still want that).
#
# Usage:  run_ld_calculation.sh -i <input_dir> -o <output_dir> -p <cpus>
################################################################################

usage() {
  echo "Usage: $0 -i <input_fasta_dir> -o <output_dir> -p <cpus>"
  exit 1
}

# Parse input arguments
while getopts "i:o:p:" opt; do
  case "$opt" in
    i) input_dir="$OPTARG";;
    o) output_dir="$OPTARG";;
    p) cpus="$OPTARG";;
    *) usage;;
  esac
done

# Check required arguments
if [ -z "${input_dir:-}" ] || [ -z "${output_dir:-}" ] || [ -z "${cpus:-}" ]; then
  usage
fi

# Ensure input directory exists
if [ ! -d "$input_dir" ]; then
  echo "ERROR: Input directory does not exist: $input_dir"
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

echo "Input directory  : $input_dir"
echo "Output directory : $output_dir"
echo "CPUs allocated   : $cpus"
echo

################################################################################
# 1. Convert Pirate's Alignment to VCF (If Needed)
################################################################################

core_alignment="$input_dir/core_alignment.fasta"  # Adjust path/name if needed
if [ ! -f "$core_alignment" ]; then
  echo "ERROR: Core alignment file not found: $core_alignment"
  exit 1
fi

echo "Checking for snp-sites..."
if ! command -v snp-sites &>/dev/null; then
  echo "WARNING: 'snp-sites' not found in PATH. Skipping snp-sites steps."
  exit 1
fi

vcf_full="$output_dir/core_snps.vcf"
if [ -f "$vcf_full" ]; then
  echo "WARNING: $vcf_full already exists, skipping snp-sites..."
else
  echo "Converting Pirate's core alignment to VCF with snp-sites..."
  snp-sites -v -o "$vcf_full" "$core_alignment"
  echo "VCF generation completed: $vcf_full"
fi

################################################################################
# 2. Randomly Subsample ~5% of Variants from the VCF
################################################################################
# This step drastically reduces the number of variants before PLINK to speed up LD
# calculations. Adjust the fraction if desired.

vcf_sub="$output_dir/core_snps_subsampled.vcf"
if [ -f "$vcf_sub" ]; then
  echo "WARNING: $vcf_sub already exists, skipping subsampling..."
else
  echo "Subsampling ~5% of variants from $vcf_full..."
  # 1) Copy header lines (start with '#')
  grep '^#' "$vcf_full" > "$vcf_sub"

  # 2) Randomly select 5% of non-header (variant) lines
  grep -v '^#' "$vcf_full" | awk 'BEGIN {srand()} rand() <= 0.05' >> "$vcf_sub"

  echo "Subsampled VCF created: $vcf_sub"
fi

################################################################################
# 3. (Optional) Further Filter the VCF with bcftools
################################################################################
# If you want additional quality-based filters, uncomment or modify below:
# bcftools view -i '%QUAL>50' -m2 -M2 -v snps -o "$output_dir/core_snps_filtered.vcf" "$vcf_sub"
# bcftools index "$output_dir/core_snps_filtered.vcf"
# Then, feed the resulting file into PLINK instead of $vcf_sub.

# For now, we just rename:
cp "$vcf_sub" "$output_dir/core_filtered.vcf"
echo "Using core_filtered.vcf for PLINK."
echo

################################################################################
# 4. Convert the Subsampled VCF to PLINK Format
################################################################################
plink_bin="./plink_linux_x86_64_20231211/plink"  # Adjust path to your PLINK binary

if [ ! -x "$plink_bin" ]; then
  echo "ERROR: PLINK binary not found or not executable at $plink_bin"
  exit 1
fi

echo "Converting core_filtered.vcf to PLINK binary format..."
"$plink_bin" --vcf "$output_dir/core_filtered.vcf" \
             --allow-extra-chr \
             --double-id \
             --make-bed \
             --out "$output_dir/plink_data"

echo "PLINK binary files created: plink_data.bed, plink_data.bim, plink_data.fam"
echo

################################################################################
# 5. Calculate LD with PLINK
################################################################################
echo "Calculating pairwise LD (R²) in PLINK..."
"$plink_bin" --bfile "$output_dir/plink_data" \
             --allow-extra-chr \
             --ld-window 200000 \
             --ld-window-kb 3000 \
             --ld-window-r2 0 \
             --r2 \
             --out "$output_dir/ld_results"

ld_file="$output_dir/ld_results.ld"
if [ -f "$ld_file" ]; then
  echo "LD results file found: $ld_file"
else
  echo "ERROR: $ld_file not produced. Check PLINK logs."
  exit 1
fi

echo

################################################################################
# 6. (Optional) Subsample LD Results for Plotting
################################################################################
# If you still want to visualize just 5% of the pairwise LD entries:

ld_sampled="$output_dir/ld_results_sampled.ld"
echo "Subsampling 5% of LD pairs for plotting..."
head -n 1 "$ld_file" > "$ld_sampled"  # Preserve header
tail -n +2 "$ld_file" | awk 'BEGIN {srand()} rand() <= 0.05' >> "$ld_sampled"

echo "Subsampled LD results: $ld_sampled"
echo

################################################################################
# 7. Cleanup Temporary Files
################################################################################
echo "Cleaning up temporary PLINK files..."
rm -f "$output_dir/plink_data."*

echo
echo "LD pipeline completed successfully!"
echo "Key outputs in $output_dir:"
echo " - core_snps_subsampled.vcf : ~5% subsample of the original Pirate VCF"
echo " - core_filtered.vcf        : VCF used for PLINK (identical to subsampled here)
echo " - plink_data.*             : PLINK binary files
echo " - ld_results.ld            : Full LD results (PLINK)
echo " - ld_results_sampled.ld    : 5% subset of LD pairs for plotting
echo
