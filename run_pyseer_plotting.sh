#!/bin/bash

# Exit if any command fails
set -e

###############################################################################
# Usage
###############################################################################
usage() {
    echo "Usage: $0 -r <ref_fna_dir> -g <gff_dir> -s <pyseer_scripts_dir> -o <output_dir> -p <prefix> [-f <specific_fasta>] [-h]"
    echo
    echo "This script will:"
    echo "  (1) Create a file named {prefix}_reference.tsv by pairing .fna/.fas/.fasta with .gff3/.gff."
    echo "  (2) Run phandango_mapper-runner.py using either:"
    echo "       (2.1) The first .fna/.fas/.fasta if no -f is provided, or"
    echo "       (2.2) The user-provided FASTA (-f) if specified."
    echo "       The output plot will have the FASTA filename appended as a suffix."
    echo "  (3) Use {prefix}_reference.tsv to run annotate_hits_pyseer-runner.py."
    echo "  (4) Summarize low/high outlier k-mers annotation results."
    echo
    echo "Required arguments:"
    echo "  -r   Directory containing .fna/.fas/.fasta files"
    echo "  -g   Directory containing .gff or .gff3 files"
    echo "  -s   Directory containing Pyseer scripts (phandango_mapper-runner.py, etc.)"
    echo "  -o   Output directory"
    echo "  -p   Prefix to name your {prefix}_reference.tsv file and some results"
    echo
    echo "Optional arguments:"
    echo "  -f   Specific FASTA file (within -r) to use for phandango_mapper-runner.py"
    echo "       (if omitted, uses the first .fna/.fas/.fasta found in -r)"
    echo "  -h   Show this help and exit"
    echo
    exit 1
}

###############################################################################
# Parse Arguments
###############################################################################
# Initialize variables
specific_fasta=""
while getopts ":r:g:s:o:p:f:h" opt; do
    case "$opt" in
        r ) ref_fna_dir="$OPTARG" ;;
        g ) gff_dir="$OPTARG" ;;
        s ) pyseer_script_dir="$OPTARG" ;;
        o ) output_dir="$OPTARG" ;;
        p ) prefix="$OPTARG" ;;
        f ) specific_fasta="$OPTARG" ;;
        h ) usage ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done
shift $((OPTIND - 1))

###############################################################################
# Validate Required Arguments
###############################################################################
if [ -z "$ref_fna_dir" ] || [ -z "$gff_dir" ] || [ -z "$pyseer_script_dir" ] || \
   [ -z "$output_dir" ] || [ -z "$prefix" ]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
result_dir="$output_dir"

###############################################################################
# Step 1: Create {prefix}_reference.tsv by pairing FASTA & GFF
###############################################################################
reference_tsv="${result_dir}/${prefix}_reference.tsv"
echo "[Info] Creating '${reference_tsv}'..."

# Optional: If you want a header row, uncomment:
# echo -e "FNA_Path\tGFF_Path\tReference" > "$reference_tsv"

shopt -s nullglob
fasta_files=("$ref_fna_dir"/*.{fna,fas,fasta})
shopt -u nullglob

if [ "${#fasta_files[@]}" -eq 0 ]; then
    echo "Error: No .fna/.fas/.fasta files found in '$ref_fna_dir'." >&2
    exit 1
fi

# We'll append to reference_tsv if you want a multi-line file:
> "$reference_tsv"  # clear or create

for fna in "${fasta_files[@]}"; do
    base_fna=$(basename "$fna")
    base_noext="${base_fna%.*}"

    # Try matching .gff or .gff3
    matched_gff=""
    for ext in gff gff3; do
        candidate="${gff_dir}/${base_noext}.${ext}"
        if [ -f "$candidate" ]; then
            matched_gff="$candidate"
            break
        fi
    done

    if [ -n "$matched_gff" ]; then
        # Append line:  FNA   GFF   ref
        echo -e "${fna}\t${matched_gff}\tref" >> "$reference_tsv"
    else
        echo "Warning: No matching .gff or .gff3 for '${fna}' found in '$gff_dir'." >&2
    fi
done

echo "[Info] Finished building '${reference_tsv}'."

###############################################################################
# Step 2: phandango_mapper-runner.py
###############################################################################
# 2.1 If no -f is given, pick the first FASTA from the array
# 2.2 If -f is given, that is the chosen FASTA
echo "[Info] Determining which FASTA to use for phandango_mapper-runner.py..."

chosen_fasta=""
if [ -z "$specific_fasta" ]; then
    # (2.1) Use the first file found
    chosen_fasta="${fasta_files[0]}"
    echo "[Info] No -f provided. Using the first FASTA: $chosen_fasta"
else
    # (2.2) Use the user-specified FASTA
    chosen_fasta="$ref_fna_dir/$specific_fasta"
    if [ ! -f "$chosen_fasta" ]; then
        echo "Error: '$chosen_fasta' not found (from -f argument)." >&2
        exit 1
    fi
    echo "[Info] Using user-specified FASTA: $chosen_fasta"
fi

# We'll extract the base name for suffix usage in the output
chosen_fasta_bn=$(basename "$chosen_fasta")
chosen_fasta_noext="${chosen_fasta_bn%.*}"

# phandango_mapper script path
phandango_mapper_runner="${pyseer_script_dir}/../phandango_mapper-runner.py"
if [ ! -f "$phandango_mapper_runner" ]; then
    echo "Error: phandango_mapper-runner.py not found at '$phandango_mapper_runner'." >&2
    exit 1
fi

# Example: Suppose we have a file named sccmec_kmers.filtered.txt and high_outliers_kmers.txt
# We'll run phandango for each, producing plots that end with the chosen FASTA's name
if [ -f "${result_dir}/${prefix}_kmers.filtered.txt" ]; then
    out_plot="${result_dir}/${prefix}_kmers_minimum_${chosen_fasta_noext}.plot"
    echo "[Info] Running phandango_mapper-runner.py on '${prefix}_kmers.filtered.txt'..."
    python "$phandango_mapper_runner" \
        "${result_dir}/${prefix}_kmers.filtered.txt" \
        "$chosen_fasta" \
        "$out_plot"
    echo "[Info] Generated: $out_plot"
fi

if [ -f "${result_dir}/high_outliers_kmers.txt" ]; then
    out_plot="${result_dir}/${prefix}_kmers_high_${chosen_fasta_noext}.plot"
    echo "[Info] Running phandango_mapper-runner.py on 'high_outliers_kmers.txt'..."
    python "$phandango_mapper_runner" \
        "${result_dir}/high_outliers_kmers.txt" \
        "$chosen_fasta" \
        "$out_plot"
    echo "[Info] Generated: $out_plot"
fi

###############################################################################
# Step 3: Use {prefix}_reference.tsv to run annotate_hits_pyseer-runner.py
###############################################################################
annotate_hits_runner="${pyseer_script_dir}/../annotate_hits_pyseer-runner.py"
if [ ! -f "$annotate_hits_runner" ]; then
    echo "Error: annotate_hits_pyseer-runner.py not found at '$annotate_hits_runner'." >&2
    exit 1
fi

# We'll annotate with the same low/high k-mers you mentioned
if [ -f "${result_dir}/filtered_kmers_normal_and_low.txt" ]; then
    low_annotated_out="${result_dir}/low_annotated_kmers_${prefix}.txt"
    echo "[Info] Annotating normal/low k-mers with $reference_tsv..."
    python "$annotate_hits_runner" \
        "${result_dir}/filtered_kmers_normal_and_low.txt" \
        "$reference_tsv" \
        "$low_annotated_out"
    echo "[Info] Created $low_annotated_out"
fi

if [ -f "${result_dir}/high_outliers_kmers.txt" ]; then
    high_annotated_out="${result_dir}/high_annotated_kmers_${prefix}.txt"
    echo "[Info] Annotating high outliers k-mers with $reference_tsv..."
    python "$annotate_hits_runner" \
        "${result_dir}/high_outliers_kmers.txt" \
        "$reference_tsv" \
        "$high_annotated_out"
    echo "[Info] Created $high_annotated_out"
fi

###############################################################################
# Step 4: Summarize low/high gene hits
###############################################################################
summarise_annotations="${pyseer_script_dir}/summarise_annotations.py"
if [ ! -f "$summarise_annotations" ]; then
    echo "Error: summarise_annotations.py not found at '$summarise_annotations'." >&2
    exit 1
fi

# Summarize low hits
if [ -f "${result_dir}/low_annotated_kmers_${prefix}.txt" ]; then
    out_low_genes="${result_dir}/low_gene_hits_${prefix}.txt"
    echo "[Info] Summarizing low annotated k-mers into $out_low_genes..."
    python "$summarise_annotations" \
        "${result_dir}/low_annotated_kmers_${prefix}.txt" \
        > "$out_low_genes"
    echo "[Info] Created $out_low_genes"
fi

# Summarize high hits
if [ -f "${result_dir}/high_annotated_kmers_${prefix}.txt" ]; then
    out_high_genes="${result_dir}/high_gene_hits_${prefix}.txt"
    echo "[Info] Summarizing high annotated k-mers into $out_high_genes..."
    python "$summarise_annotations" \
        "${result_dir}/high_annotated_kmers_${prefix}.txt" \
        > "$out_high_genes"
    echo "[Info] Created $out_high_genes"
fi

echo "[Success] Pipeline completed for prefix '$prefix'."
