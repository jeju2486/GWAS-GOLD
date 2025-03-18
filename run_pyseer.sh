#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -t <tree_file> -p <phenotype_file> -i <unitig_output_dir> -o <output_dir> -P <prefix> -T <threads> -s <pyseer_scripts_dir>"
    exit 1
}

# Parse command-line arguments
while getopts ":t:p:i:o:P:T:s:" opt; do
  case $opt in
    t) treefile="$OPTARG" ;;
    p) phenotype="$OPTARG" ;;
    i) unitig_result="$OPTARG" ;;
    o) result_dir="$OPTARG" ;;
    P) prefix="$OPTARG" ;;
    T) threads="$OPTARG" ;;
    s) pyseer_scripts_dir="$OPTARG" ;;
    \?)
      echo "Invalid option -$OPTARG" >&2
      usage
      ;;
  esac
done

# Check that all required arguments are provided
if [ -z "$treefile" ] || [ -z "$phenotype" ] || [ -z "$unitig_result" ] || \
   [ -z "$result_dir" ] || [ -z "$prefix" ] || [ -z "$threads" ] || [ -z "$pyseer_scripts_dir" ]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$result_dir"

# Check if the main Pyseer output file already exists
if [ ! -f "$result_dir/${prefix}_kmers.txt" ]; then
    echo "[INFO] Generating phylogeny-based similarity matrix..."
    python "$pyseer_scripts_dir"/phylogeny_distance.py \
        --lmm "$treefile" \
        > "$result_dir"/"${prefix}_phylogeny.tsv" 2>/dev/null

    echo "[INFO] Running Pyseer (this may take some time)..."
    pyseer \
        --lmm \
        --phenotypes "$phenotype" \
        --kmers "$unitig_result"/survived_unitigs.pyseer.gz \
        --similarity "$result_dir"/"${prefix}_phylogeny.tsv" \
        --min-af 0.05 \
        --output-patterns "$result_dir"/kmer_patterns.txt \
        --cpu "$threads" \
        > "$result_dir"/"${prefix}_kmers.txt" 2>/dev/null

    echo "[INFO] Pyseer analysis complete."
else
    echo "[INFO] Pyseer output file ${prefix}_kmers.txt already exists. Skipping Pyseer analysis."
fi

echo "[DONE] Results are in: $result_dir/${prefix}_kmers.txt"
