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
    echo "  (2) Run phandango_mapper-runner.py on {prefix}_kmers.txt (if present),"
    echo "      using either the first .fna/.fas/.fasta if no -f is provided, or a user-specified FASTA."
    echo "  (3) Generate a Q-Q plot using {prefix}_kmers.txt."
    echo "  (4) Annotate and summarize only significant_kmers.txt (renamed outputs)."
    echo "      *All low/high outlier steps are removed.*"
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

# Clear or create the reference file
> "$reference_tsv"

shopt -s nullglob
fasta_files=("$ref_fna_dir"/*.{fna,fas,fasta})
shopt -u nullglob

if [ "${#fasta_files[@]}" -eq 0 ]; then
    echo "Error: No .fna/.fas/.fasta files found in '$ref_fna_dir'." >&2
    exit 1
fi

# Pair each FASTA with its matching .gff or .gff3
for fna in "${fasta_files[@]}"; do
    base_fna=$(basename "$fna")
    base_noext="${base_fna%.*}"

    matched_gff=""
    for ext in gff gff3; do
        candidate="${gff_dir}/${base_noext}.${ext}"
        if [ -f "$candidate" ]; then
            matched_gff="$candidate"
            break
        fi
    done

    if [ -n "$matched_gff" ]; then
        echo -e "${fna}\t${matched_gff}\tref" >> "$reference_tsv"
    else
        echo "Warning: No matching .gff or .gff3 for '${fna}' found in '$gff_dir'." >&2
    fi
done

echo "[Info] Finished building '${reference_tsv}'."

###############################################################################
# Step 2: phandango_mapper-runner.py on {prefix}_kmers.txt
###############################################################################
# Decide which FASTA to use
echo "[Info] Determining which FASTA to use for phandango_mapper-runner.py..."
chosen_fasta=""
if [ -z "$specific_fasta" ]; then
    chosen_fasta="${fasta_files[0]}"
    echo "[Info] No -f provided. Using the first FASTA: $chosen_fasta"
else
    chosen_fasta="$ref_fna_dir/$specific_fasta"
    if [ ! -f "$chosen_fasta" ]; then
        echo "Error: '$chosen_fasta' not found (from -f argument)." >&2
        exit 1
    fi
    echo "[Info] Using user-specified FASTA: $chosen_fasta"
fi

# phandango_mapper script path
phandango_mapper_runner="${pyseer_script_dir}/../phandango_mapper-runner.py"
if [ ! -f "$phandango_mapper_runner" ]; then
    echo "Error: phandango_mapper-runner.py not found at '$phandango_mapper_runner'." >&2
    exit 1
fi

if [ -f "${result_dir}/${prefix}_kmers.txt" ]; then
    chosen_fasta_bn=$(basename "$chosen_fasta")
    chosen_fasta_noext="${chosen_fasta_bn%.*}"
    out_plot="${result_dir}/${prefix}_kmers_${chosen_fasta_noext}.plot"

    echo "[Info] Running phandango_mapper-runner.py on '${prefix}_kmers.txt'..."
    python "$phandango_mapper_runner" \
        "${result_dir}/${prefix}_kmers.txt" \
        "$chosen_fasta" \
        "$out_plot"
    echo "[Info] Generated: $out_plot"
else
    echo "[Info] '${prefix}_kmers.txt' not found; skipping phandango plotting."
fi

###############################################################################
# Step 2.5: Generate the Q-Q plot (using {prefix}_kmers.txt)
###############################################################################
qq_plot="${pyseer_script_dir}/qq_plot.py"

if [ -f "${result_dir}/${prefix}_kmers.txt" ]; then
    echo "[Info] Generating QQ-plot from ${prefix}_kmers.txt..."
    python "$qq_plot" "${result_dir}/${prefix}_kmers.txt"
    mv ./qq_plot.png "${result_dir}"
    echo "[Info] QQ-plot saved to ${result_dir}/qq_plot.png"
else
    echo "[Info] Cannot generate QQ-plot; ${prefix}_kmers.txt does not exist."
fi

###############################################################################
# Step 3: Annotate ONLY significant_kmers.txt (no more 'original_' prefix)
###############################################################################
annotate_hits_runner="${pyseer_script_dir}/../annotate_hits_pyseer-runner.py"
if [ ! -f "$annotate_hits_runner" ]; then
    echo "Error: annotate_hits_pyseer-runner.py not found at '$annotate_hits_runner'." >&2
    exit 1
fi

if [ -f "${result_dir}/significant_kmers.txt" ]; then
    annotated_out="${result_dir}/annotated_kmers_${prefix}.txt"
    echo "[Info] Annotating significant k-mers with $reference_tsv..."
    python "$annotate_hits_runner" \
        "${result_dir}/significant_kmers.txt" \
        "$reference_tsv" \
        "$annotated_out"
    echo "[Info] Created $annotated_out"
fi

###############################################################################
# Step 4: Summarize only the annotated significant_kmers (no low/high blocks)
###############################################################################
summarise_annotations="${pyseer_script_dir}/summarise_annotations.py"
if [ ! -f "$summarise_annotations" ]; then
    echo "Error: summarise_annotations.py not found at '$summarise_annotations'." >&2
    exit 1
fi

if [ -f "${result_dir}/annotated_kmers_${prefix}.txt" ]; then
    out_genes="${result_dir}/gene_hits_${prefix}.txt"
    echo "[Info] Summarizing annotated significant k-mers into $out_genes..."
    python "$summarise_annotations" --nearby \
        "${result_dir}/annotated_kmers_${prefix}.txt" \
        > "$out_genes"
    echo "[Info] Created $out_genes"
fi

###############################################################################
# Done
###############################################################################
echo "[Success] Pipeline completed for prefix '$prefix'."
