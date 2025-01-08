################################################################################
# 1) PARSE COMMAND-LINE ARGUMENTS
################################################################################

usage() {
    echo "Usage: $0 -i <input_dir> -d <ld_length> -o <output_dir> -t <thread_num> [-q <query_sequence> | -s <summary_tsv>]"
    echo
    echo "Options:"
    echo "  -q    Query sequence to map (FASTA file)"
    echo "  -s    Summary TSV file with pre-computed results"
    echo "  -i    Directory containing reference FASTA files"
    echo "  -d    Length for LD region"
    echo "  -o    Directory to store outputs"
    echo "  -t    Number of threads"
    echo "  -h    Show this help message"
    exit 1
}

# Initialize variables
summary_tsv=""
query_sequence=""

# Parse arguments
while getopts ":q:s:i:d:o:t:h" flag
do
    case "${flag}" in
        q) query_sequence=${OPTARG};;    # Query sequence to map
        s) summary_tsv=${OPTARG};;       # TSV with pre-computed results
        i) input_dir=${OPTARG};;         # Directory containing reference FASTA
        d) ld_length=${OPTARG};;         # Length for LD region
        o) output_dir=${OPTARG};;        # Where to store outputs
        t) thread_num=${OPTARG};;        # Number of threads
        h) usage;;                       # Show help message
        \?) echo "Invalid option: -$OPTARG" >&2; usage;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    esac
done

# Validate required arguments
if [[ -z "${input_dir}" || -z "${ld_length}" || -z "${output_dir}" || -z "${thread_num}" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Ensure that either -q or -s is provided, but not both
if [[ -n "${query_sequence}" && -n "${summary_tsv}" ]]; then
    echo "Error: Specify either -q (query sequence) or -s (summary TSV), not both."
    usage
elif [[ -z "${query_sequence}" && -z "${summary_tsv}" ]]; then
    echo "Error: You must specify either -q (query sequence) or -s (summary TSV)."
    usage
fi

# Check if provided files/directories exist
if [[ -n "${query_sequence}" && ! -f "${query_sequence}" ]]; then
    echo "Error: Query sequence file '${query_sequence}' does not exist."
    exit 1
fi

if [[ -n "${summary_tsv}" && ! -f "${summary_tsv}" ]]; then
    echo "Error: Summary TSV file '${summary_tsv}' does not exist."
    exit 1
fi

if [[ ! -d "${input_dir}" ]]; then
    echo "Error: Input directory '${input_dir}' does not exist."
    exit 1
fi

################################################################################
# 2) CREATE NECESSARY DIRECTORIES
################################################################################
mkdir -p "${output_dir}/sam" \
         "${output_dir}/bed" \
         "${output_dir}/ld_ref" \
         "${output_dir}/temp" \
         "${output_dir}/unitig_output"

# Remove any old files in ld_ref (if any exist)
rm -f "${output_dir}/ld_ref/"*

################################################################################
# 3) OPTIONAL: PARSE THE SUMMARY TSV INTO ARRAYS/DICTIONARIES (IF PROVIDED)
################################################################################
declare -A ref2contig
declare -A ref2start
declare -A ref2end

if [[ -n "${summary_tsv}" ]]; then
    echo "Reading summary TSV: ${summary_tsv}"
    # Expected TSV format: Reference Contig Start0 End
    # Example:
    # GCF_000144955.fas NC_017338.2 33707 57871
    {
        read -r header  # skip header
        while read -r reference contig start0 endpos; do
            # Store only if reference is not empty and coordinates are valid
            if [[ -n "${reference}" && "${contig}" != "NA" && "${start0}" != "NA" && "${endpos}" != "NA" ]]; then
                ref2contig["${reference}"]="${contig}"
                ref2start["${reference}"]="${start0}"
                ref2end["${reference}"]="${endpos}"
            fi
        done
    } < "${summary_tsv}"
fi

################################################################################
# 4) ALIGN QUERY SEQUENCE OR USE TSV DATA TO CREATE BED FILES
################################################################################

for reference_genome in "${input_dir}"/*.{fas,fasta,fna}; do

    # Skip if no matching file
    if [[ ! -e $reference_genome ]]; then
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"

    # -----------------------------------------------------------------------------
    # A) If summary TSV is provided and contains data for this reference
    # -----------------------------------------------------------------------------
    if [[ -n "${summary_tsv}" && -n "${ref2contig["${reference_name}"]}" ]]; then
        # Pre-computed data exists
        contig="${ref2contig["${reference_name}"]}"
        start0="${ref2start["${reference_name}"]}"
        endpos="${ref2end["${reference_name}"]}"

        echo "Using TSV data for reference: ${reference_name}"
        echo "Contig = ${contig}, Start0 = ${start0}, End = ${endpos}"

        # Create a pseudo-SAM file if downstream steps require it
        sam_file="${output_dir}/sam/${reference_base}.sam"
        if [ ! -e "${sam_file}" ]; then
            echo -e "@HD\tVN:1.6\tSO:unsorted" > "${sam_file}"
            echo -e "@SQ\tSN:${contig}\tLN:999999999" >> "${sam_file}"
            # Placeholder alignment entry
            echo -e "READ_1\t0\t${contig}\t${start0}\t60\t*\t*\t0\t0\t*\t*\tAS:i:0\tXS:i:0" \
                >> "${sam_file}"
        else
            echo "SAM file '${sam_file}' already exists. Skipping SAM creation."
        fi

        # Create BED file from TSV coordinates
        bed_file="${output_dir}/bed/${reference_base}.bed"
        if [ ! -s "${bed_file}" ]; then
            # Assuming '+' strand; adjust if strand information is available
            echo -e "${contig}\t${start0}\t${endpos}\t${reference_base}\t+" \
                > "${bed_file}"
        else
            echo "BED file '${bed_file}' already exists. Skipping BED creation."
        fi

    else
        # -----------------------------------------------------------------------------
        # B) Fallback: No pre-computed data or summary TSV not provided => run minimap2
        # -----------------------------------------------------------------------------
        sam_file="${output_dir}/sam/${reference_base}.sam"
        bed_file="${output_dir}/bed/${reference_base}.bed"

        if [ -e "${sam_file}" ]; then
            echo "SAM file '${sam_file}' already exists. Skipping minimap2."
        else
            if [[ -n "${query_sequence}" ]]; then
                echo "Processing reference (minimap2): ${reference_base}"
                # Run Minimap2 and save the SAM, suppress standard error
                minimap2 -a -t "${thread_num}" "${reference_genome}" "${query_sequence}" \
                    2>/dev/null \
                    > "${sam_file}"
            else
                echo "Error: No query sequence provided and no summary data available for '${reference_name}'."
                echo "Skipping this reference."
                continue
            fi
        fi

        # Convert SAM to BED if BED doesn't already exist
        if [ -e "${sam_file}" ] && [ ! -s "${bed_file}" ]; then
            bedtools bamtobed -i "${sam_file}" \
                > "${bed_file}"
        fi
    fi

    # Create genome length information file (for slop steps) if not already done
    genome_size_file="${output_dir}/bed/${reference_base}_genome_size.txt"
    if [ ! -s "${genome_size_file}" ]; then
        awk '/^>/{                                \
                if (seqname)                      \
                    print seqname "\t" length(seq);\
                seqname=$1;                       \
                seq="";                           \
                next                              \
             }                                    \
             {seq = seq $0}                      \
             END {print seqname "\t" length(seq)}' "${reference_genome}" \
          | sed 's/>//' \
          > "${genome_size_file}"
    fi

done

################################################################################
# 5) CREATE EXTENDED BED REGIONS FOR UPSTREAM/DOWNSTREAM LD, AND GET FASTA
################################################################################
for reference_genome in "${input_dir}"/*.{fas,fasta,fna}; do

    # Skip if no matching file
    if [[ ! -e $reference_genome ]]; then
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"

    echo "Processing extended LD for: ${reference_base}"

    bed_file="${output_dir}/bed/${reference_base}.bed"
    genome_size_file="${output_dir}/bed/${reference_base}_genome_size.txt"

    # Skip if BED file does not exist or is empty
    if [ ! -s "${bed_file}" ]; then
        echo "No gene found (BED empty) for '${reference_base}', Skipping."
        continue
    fi

    # Initialize run counter
    run_num=1

    # Read each line in BED
    while IFS=$'\t' read -r chr start end gene_name strand; do
        # Calculate gene length
        query_length=$((end - start))
        echo "  Creating LD for gene: ${gene_name}, #${run_num}, length: ${query_length}"

        # Temporary BED for single gene
        temp_bed="${output_dir}/temp/${reference_base}_${run_num}.bed"
        echo -e "${chr}\t${start}\t${end}\t${gene_name}\t${strand}" \
            > "${temp_bed}"

        # Extend regions by ld_length on each side
        bedtools slop \
            -i "${temp_bed}" \
            -g "${genome_size_file}" \
            -b "${ld_length}" \
            >> "${output_dir}/ld_ref/${reference_base}_elongated.bed"

        run_num=$((run_num + 1))
    done < "${bed_file}"

    # Extract fasta for the elongated (LD) regions
    ld_fasta="${output_dir}/ld_ref/${reference_base}_query_gene.fas"
    bedtools getfasta \
        -fi "${reference_genome}" \
        -bed "${output_dir}/ld_ref/${reference_base}_elongated.bed" \
        -fo "${ld_fasta}"
done

################################################################################
# 6) CLEAN UP AND PREPARE LISTS FOR UNITIG-CALLER
################################################################################
echo "Removing temporary files..."
rm -rf "${output_dir}/temp"
rm -f "${input_dir}"/*.fai \
      "${input_dir}"/*.amb \
      "${input_dir}"/*.ann \
      "${input_dir}"/*.bwt \
      "${input_dir}"/*.pac \
      "${input_dir}"/*.sa

# Find .fas files in ld_ref and write to target_input.txt
find "${output_dir}/ld_ref/" -maxdepth 1 -type f -name "*.fas" -print \
    > "${output_dir}/target_input.txt"

# Find FASTA files in input_dir and write to genome_input.txt
find "${input_dir}" -maxdepth 1 -type f \( -name "*.fas" -o -name "*.fasta" -o -name "*.fna" \) -print \
    > "${output_dir}/genome_input.txt"

################################################################################
# 7) RUN UNITIG-CALLER
################################################################################
unitig_output_pyseer="${output_dir}/unitig_output/unitig.out.pyseer"

if [ -f "${unitig_output_pyseer}" ]; then
    echo "unitig output '${unitig_output_pyseer}' already exists. Skipping unitig-caller."
else
    echo "Running unitig-caller for genomes..."
    unitig-caller --call \
        --refs "${output_dir}/genome_input.txt" \
        --out  "${output_dir}/unitig_output/unitig.out" \
        --kmer 31 \
        --pyseer \
        --threads "${thread_num}"

    gzip "${output_dir}/unitig_output/extracted.unitig.out.pyseer"
    echo "unitig-caller finished"
fi

################################################################################
# 8) CONVERT UNITIGS TO FASTA
################################################################################
unitig_fasta="${output_dir}/unitig_output/unitig.out.fasta"

awk -F'|' '{print ">" NR "\n" $1}' "${unitig_output_pyseer}" \
    > "${unitig_fasta}"

################################################################################
# 9) FILTER UNITIGS AGAINST COMBINED REFERENCE SEQUENCES USING BWA
################################################################################
bwa_output_dir="bwa_output"
mkdir -p "${output_dir}/unitig_output/${bwa_output_dir}"

# Combine relevant reference sequences
combined_fasta="${output_dir}/unitig_output/combined_sequences.fasta"
cat "${output_dir}/ld_ref/"*.fas \
    > "${combined_fasta}"

# Index combined sequences
bwa index "${combined_fasta}" 2>/dev/null

# Align unitigs to combined sequences, suppress BWA printing
bwa mem -k 13 -t "${thread_num}" \
    "${combined_fasta}" \
    "${unitig_fasta}" \
    2>/dev/null \
    > "${output_dir}/unitig_output/${bwa_output_dir}/output.sam"

# Convert SAM to BED
bedtools bamtobed \
    -i "${output_dir}/unitig_output/${bwa_output_dir}/output.sam" \
    > "${output_dir}/unitig_output/${bwa_output_dir}/output.bed"

# Gather IDs that actually aligned
awk -F'\t' '!/^@/ && $3 != "*" {print $1}' \
    "${output_dir}/unitig_output/${bwa_output_dir}/output.sam" \
    | sort -n | uniq \
    > "${output_dir}/unitig_output/ids_to_extract.txt"

# Filter out matched unitigs
awk 'NR==FNR {a[$1]; next} /^>/ {header=$1; sub(/^>/, "", header); keep = !(header in a)} keep' \
    "${output_dir}/unitig_output/ids_to_extract.txt" \
    "${unitig_fasta}" \
    > "${output_dir}/unitig_output/unitig.filtered.fasta"

# Extract only the sequences (not headers)
grep -v '^>' "${output_dir}/unitig_output/unitig.filtered.fasta" \
    > "${output_dir}/unitig_output/survived_sequences.txt"

# Remove empty lines
sed -i '/^$/d' "${output_dir}/unitig_output/survived_sequences.txt"

################################################################################
# 10) GREP CORRESPONDING LINES FROM PYSEER UNITIG OUTPUT
################################################################################
survived_sequences="${output_dir}/unitig_output/survived_sequences.txt"
survived_unitigs="${output_dir}/unitig_output/survived_unitigs.pyseer"

grep -F -f "${survived_sequences}" "${unitig_output_pyseer}" \
    > "${survived_unitigs}"

# Gzip the final set
survived_unitigs_gz="${survived_unitigs}.gz"
if [ -e "${survived_unitigs_gz}" ]; then
    rm "${survived_unitigs_gz}"
fi
gzip "${survived_unitigs}"

################################################################################
# 11) FINAL CLEANUP
################################################################################
echo "Removing temporary unitig files..."
rm -f "${unitig_fasta}" \
      "${output_dir}/unitig_output/unitig.filtered.fasta" \
      "${survived_sequences}" \
      "${output_dir}/unitig_output/ids_to_extract.txt"

echo "All tasks completed successfully."
echo "Ready to run pyseer or any further analyses."
