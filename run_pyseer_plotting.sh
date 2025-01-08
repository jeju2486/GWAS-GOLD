#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

################################################################################
# Logging (Optional)
################################################################################
# Uncomment the following lines to enable logging to a file named run_pyseer_plotting.log
# located in the output directory.
# Make sure to place these lines after the output directory is created.

# exec > >(tee -i "$output_dir/run_pyseer_plotting.log")
# exec 2>&1

################################################################################
# Usage Function
################################################################################
usage() {
    echo "Usage: $0 -r <ref_fna_dir> -g <gff_dir> [-e <existing_ref_file>] -s <pyseer_scripts_dir> -o <output_dir> -p <prefix>"
    echo
    echo "This script will:"
    echo "  1) Generate a reference file (unless one is provided) by pairing .fna/.fas/.fasta and .gff3/.gff files."
    echo "  2) Index the reference sequence with BWA."
    echo "  3) Run phandango_mapper-runner.py on specified k-mer files."
    echo "  4) Run annotate_hits_pyseer-runner.py on k-mer files using the reference."
    echo "  5) Summarize annotation hits."
    echo
    echo "Required arguments:"
    echo "  -r   Directory containing .fna, .fas, or .fasta reference files (e.g., /path/to/filtered)"
    echo "  -g   Directory containing .gff3 or .gff files (Prokka outputs) (e.g., /path/to/prokka_gffs)"
    echo "  -s   Pyseer scripts directory (e.g., /path/to/pyseer_codes/scripts)"
    echo "  -o   Output directory (results will be saved here)"
    echo "  -p   Prefix for naming output files (e.g., GCF_000144955)"
    echo
    echo "Optional arguments:"
    echo "  -e   Existing reference file (e.g., reference_sepi_add.tsv or existing_reference.fasta). If provided, the script "
    echo "       will skip generating a new reference file and use this as the reference sequence."
    echo "  -h   Show this help message and exit"
    echo
    exit 1
}

################################################################################
# Parse Arguments
################################################################################
# Initialize variables
while getopts ":r:g:e:s:o:p:h" opt; do
    case ${opt} in
        r )
            ref_fna_dir="$OPTARG"
            ;;
        g )
            gff_dir="$OPTARG"
            ;;
        e )
            existing_ref_file="$OPTARG"
            ;;
        s )
            pyseer_script_dir="$OPTARG"
            ;;
        o )
            output_dir="$OPTARG"
            ;;
        p )
            prefix="$OPTARG"
            ;;
        h )
            usage
            ;;
        \? )
            echo "Invalid Option: -$OPTARG" >&2
            usage
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Shift off the options and optional --
shift $((OPTIND -1))

################################################################################
# Validate Required Arguments
################################################################################
if [ -z "$ref_fna_dir" ] || [ -z "$gff_dir" ] || [ -z "$pyseer_script_dir" ] || \
   [ -z "$output_dir" ] || [ -z "$prefix" ]; then
    echo "Error: Missing one or more required arguments." >&2
    usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
result_dir="$output_dir"

################################################################################
# Step 1: Generate the Reference File (if Needed)
################################################################################
if [ -z "$existing_ref_file" ]; then
    echo "[Info] No existing reference file provided. Generating a combined reference file..."

    # Define the output reference file path
    combined_ref_file="${result_dir}/reference.tsv"

    # Loop through each .fna, .fas, or .fasta file in the reference directory
    shopt -s nullglob
    fasta_files=("$ref_fna_dir"/*.{fna,fas,fasta})
    shopt -u nullglob

    if [ "${#fasta_files[@]}" -eq 0 ]; then
        echo "Error: No FASTA files with extensions .fna, .fas, or .fasta found in '$ref_fna_dir'." >&2
        exit 1
    fi

    for fna in "${fasta_files[@]}"; do
        # Extract the base filename without extension
        base=$(basename "$fna")
        base_no_ext="${base%.*}"

        # Search for corresponding .gff3 or .gff file
        gff_found=0
        for gff_ext in gff3 gff; do
            gff_file="${gff_dir}/${base_no_ext}.${gff_ext}"
            if [[ -f "$gff_file" ]]; then
                # Append the paths and "ref" to the output file, separated by tabs
                echo -e "${fna}\t${gff_file}\tref" >> "$combined_ref_file"
                gff_found=1
                break
            fi
        done

        if [[ $gff_found -eq 0 ]]; then
            echo "Warning: Corresponding .gff3 or .gff file for '${fna}' not found." >&2
        fi
    done

    echo "[Info] Reference file '${combined_ref_file}' has been created successfully."

    # Determine the reference_seq to use for phandango_mapper-runner.py
    # Use the first FASTA file as the reference sequence
    reference_seq="${fasta_files[0]}"
    echo "[Info] Selected reference sequence: '$reference_seq'"

else
    # Use the existing reference file provided by the user
    echo "[Info] Using existing reference file: '${existing_ref_file}'"
    reference_seq="$existing_ref_file"

    # Validate the existence of the provided reference file
    if [[ ! -f "$reference_seq" ]]; then
        echo "Error: Provided reference file '$reference_seq' does not exist." >&2
        exit 1
    fi

    # If the existing reference file is a TSV (e.g., reference_sepi_add.tsv), select the first FASTA entry
    # Assuming the third column 'Reference' contains the identifier 'ref' for the reference
    # Modify this part as per your actual reference file structure if needed

    # Check if the existing_ref_file is a TSV (contains tab-separated values)
    if [[ "$reference_seq" == *.tsv ]]; then
        echo "[Info] Existing reference file is a TSV. Extracting the first FASTA file as the reference sequence."

        # Extract the first FASTA path from the TSV
        first_fasta=$(awk -F'\t' '$3 == "ref" {print $1; exit}' "$reference_seq")

        if [[ -z "$first_fasta" ]]; then
            echo "Error: No reference FASTA entry found in '$reference_seq'." >&2
            exit 1
        fi

        reference_seq="$first_fasta"
        echo "[Info] Selected reference sequence from TSV: '$reference_seq'"
    fi
fi

################################################################################
# Step 1.5: Index the Reference Sequence with BWA
################################################################################
echo "[Info] Indexing reference sequence with BWA..."
bwa index "$reference_seq"
echo "[Info] BWA indexing completed."

################################################################################
# Step 2: Define Paths to Pyseer Scripts
################################################################################
phandango_mapper_runner="${pyseer_script_dir}/../phandango_mapper-runner.py"
annotate_hits_runner="${pyseer_script_dir}/../annotate_hits_pyseer-runner.py"
summarise_annotations="${pyseer_script_dir}/summarise_annotations.py"

# Validate that the scripts exist
for script in "$phandango_mapper_runner" "$annotate_hits_runner" "$summarise_annotations"; do
    if [[ ! -f "$script" ]]; then
        echo "Error: Required script '$script' not found." >&2
        exit 1
    fi
done

################################################################################
# Step 3: Check for Required Input K-mer Files
################################################################################
required_kmer_files=(
    "${result_dir}/ecoli_kmers.filtered.txt"
    "${result_dir}/high_outliers_kmers.txt"
    "${result_dir}/filtered_kmers_normal_and_low.txt"
)

missing_files=()

for file in "${required_kmer_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "Warning: Expected k-mer file '$file' not found." >&2
        missing_files+=("$file")
    fi
done

if [[ "${#missing_files[@]}" -eq "${#required_kmer_files[@]}" ]]; then
    echo "Error: All required k-mer files are missing. Exiting." >&2
    exit 1
elif [[ "${#missing_files[@]}" -gt 0 ]]; then
    echo "Proceeding with available k-mer files."
fi

################################################################################
# Step 4: Run phandango_mapper-runner.py for K-mer Files
################################################################################
# Function to run phandango_mapper-runner.py
run_phandango_mapper() {
    local kmer_file="$1"
    local reference="$2"
    local output_plot="$3"

    echo "[Info] Running phandango_mapper-runner.py for '$kmer_file'..."
    python "$phandango_mapper_runner" \
        "$kmer_file" \
        "$reference" \
        "$output_plot"
    echo "[Info] Generated plot: '$output_plot'"
}

# Run for prefix_kmers.filtered.txt
if [[ -f "${result_dir}/ecoli_kmers.filtered.txt" ]]; then
    run_phandango_mapper \
        "${result_dir}/ecoli_kmers.filtered.txt" \
        "$reference_seq" \
        "${result_dir}/ecoli_kmers_minimum_${prefix}.plot"
fi

# Run for high_outliers_kmers.txt
if [[ -f "${result_dir}/high_outliers_kmers.txt" ]]; then
    run_phandango_mapper \
        "${result_dir}/high_outliers_kmers.txt" \
        "$reference_seq" \
        "${result_dir}/ecoli_kmers_high_${prefix}.plot"
fi

################################################################################
# Step 5: Run annotate_hits_pyseer-runner.py for K-mer Files
################################################################################
# Determine the reference file for annotation
if [[ -z "$existing_ref_file" ]]; then
    annotate_ref_file="$combined_ref_file"
else
    annotate_ref_file="$existing_ref_file"
fi

# Validate that annotate_ref_file exists
if [[ ! -f "$annotate_ref_file" ]]; then
    echo "Error: Reference file for annotation '$annotate_ref_file' does not exist." >&2
    exit 1
fi

# Function to run annotate_hits_pyseer-runner.py
run_annotate_hits() {
    local kmer_file="$1"
    local reference_file="$2"
    local output_annotated_kmers="$3"

    echo "[Info] Annotating k-mers in '$kmer_file'..."
    python "$annotate_hits_runner" \
        "$kmer_file" \
        "$reference_file" \
        "$output_annotated_kmers"
    echo "[Info] Generated annotated k-mers file: '$output_annotated_kmers'"
}

# Annotate filtered_kmers_normal_and_low.txt
if [[ -f "${result_dir}/filtered_kmers_normal_and_low.txt" ]]; then
    run_annotate_hits \
        "${result_dir}/filtered_kmers_normal_and_low.txt" \
        "$annotate_ref_file" \
        "${result_dir}/low_annotated_kmers_${prefix}.txt"
fi

# Annotate high_outliers_kmers.txt
if [[ -f "${result_dir}/high_outliers_kmers.txt" ]]; then
    run_annotate_hits \
        "${result_dir}/high_outliers_kmers.txt" \
        "$annotate_ref_file" \
        "${result_dir}/high_annotated_kmers_${prefix}.txt"
fi

################################################################################
# Step 6: Summarize Annotations
################################################################################
# Function to summarize annotations
summarize_annotations() {
    local annotated_kmers_file="$1"
    local output_gene_hits="$2"

    echo "[Info] Summarizing annotations in '$annotated_kmers_file'..."
    python "$summarise_annotations" \
        "$annotated_kmers_file" \
        > "$output_gene_hits"
    echo "[Info] Generated gene hits file: '$output_gene_hits'"
}

# Summarize low_annotated_kmers
if [[ -f "${result_dir}/low_annotated_kmers_${prefix}.txt" ]]; then
    summarize_annotations \
        "${result_dir}/low_annotated_kmers_${prefix}.txt" \
        "${result_dir}/low_gene_hits_${prefix}.txt"
fi

# Summarize high_annotated_kmers
if [[ -f "${result_dir}/high_annotated_kmers_${prefix}.txt" ]]; then
    summarize_annotations \
        "${result_dir}/high_annotated_kmers_${prefix}.txt" \
        "${result_dir}/high_gene_hits_${prefix}.txt"
fi

################################################################################
# Completion Message
################################################################################
echo "[Success] Pyseer plotting pipeline completed successfully."
