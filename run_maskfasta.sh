#!/usr/bin/env bash

################################################################################
# 0) SETUP DEBUG FUNCTION
################################################################################

debug=0  # default: off

# Function to print debug messages only if debug=1
debug_msg() {
    if [[ "${debug}" -eq 1 ]]; then
        echo "$@"
    fi
}

################################################################################
# 1) PARSE COMMAND-LINE ARGUMENTS
################################################################################

usage() {
    echo "Usage: $0 -i <input_dir> -d <ld_length> -o <output_dir> -t <thread_num> [-s <summary_tsv>] -g <gene_fasta> [-x]"
    echo
    echo "Options:"
    echo "  -s    Summary TSV file with pre-computed results"
    echo "  -i    Directory containing reference FASTA files"
    echo "  -d    Length for LD region"
    echo "  -o    Directory to store outputs"
    echo "  -t    Number of threads"
    echo "  -g    Gene FASTA file (e.g., containing lukF-PV/lukS-PV), used to build ABRicate DB"
    echo "  -x    Enable debug (verbose) messages"
    echo "  -h    Show this help message"
    exit 1
}

# Initialize variables
summary_tsv=""
input_dir=""
ld_length=""
output_dir=""
thread_num=""
gene_fasta=""

# Parse arguments (note the added -x for debug)
while getopts ":s:i:d:o:t:g:hx" flag; do
    case "${flag}" in
        s) summary_tsv=${OPTARG};;   # TSV with pre-computed results
        i) input_dir=${OPTARG};;     # Directory containing reference FASTA
        d) ld_length=${OPTARG};;     # Length for LD region
        o) output_dir=${OPTARG};;    # Where to store outputs
        t) thread_num=${OPTARG};;    # Number of threads
        g) gene_fasta=${OPTARG};;    # Gene FASTA for ABRicate
        x) debug=1;;                 # Debug/verbose flag
        h) usage;;                   # Show help message
        \?) echo "Error: Invalid option -$OPTARG" >&2; usage;;
        :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage;;
    esac
done

# Validate required arguments
if [[ -z "${input_dir}" || -z "${ld_length}" || -z "${output_dir}" || -z "${thread_num}" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# If no summary TSV is given, we require a gene FASTA
if [[ -z "${summary_tsv}" && -z "${gene_fasta}" ]]; then
    echo "Error: You must provide either a summary TSV (-s) or a gene FASTA (-g) for ABRicate."
    usage
fi

# Check if provided files/directories exist
if [[ -n "${summary_tsv}" && ! -f "${summary_tsv}" ]]; then
    echo "Error: Summary TSV file '${summary_tsv}' does not exist."
    exit 1
fi

if [[ ! -d "${input_dir}" ]]; then
    echo "Error: Input directory '${input_dir}' does not exist."
    exit 1
fi

if [[ -n "${gene_fasta}" && ! -f "${gene_fasta}" ]]; then
    echo "Error: Gene FASTA file '${gene_fasta}' does not exist."
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
# 3) OPTIONAL: PARSE THE SUMMARY TSV INTO ARRAYS (IF PROVIDED)
################################################################################
declare -A ref2contig
declare -A ref2start
declare -A ref2end

if [[ -n "${summary_tsv}" ]]; then
    debug_msg "Reading summary TSV: ${summary_tsv}"
    # Expected TSV format: ReferenceFasta  Contig  Start0  End
    {
        read -r header # skip header
        while read -r reference gene_presence contig start0 endpos; do
            if [[ -n "${reference}" && "${contig}" != "NA" && "${start0}" != "NA" && "${endpos}" != "NA" ]]; then
                ref2contig["${reference}"]="${contig}"
                ref2start["${reference}"]="${start0}"
                ref2end["${reference}"]="${endpos}"
            fi
        done
    } < "${summary_tsv}"
fi

################################################################################
# 4) BUILD ABRICATE DATABASE (IF NO SUMMARY TSV PROVIDED)
################################################################################

abricate_db_dir="${output_dir}/abricate_db"
abricate_db_name="my_abricate_db"

if [[ -z "${summary_tsv}" ]]; then
    debug_msg "Making ABRicate DB from '${gene_fasta}'..."

    # 1. Create the subdirectory named exactly how you'll refer to it with --db
    mkdir -p "${abricate_db_dir}/${abricate_db_name}"

    # 2. Copy your gene FASTA to 'sequences' inside that directory
    cp "${gene_fasta}" "${abricate_db_dir}/${abricate_db_name}/sequences"

    # 3. Build the BLAST database from 'sequences'
    pushd "${abricate_db_dir}/${abricate_db_name}" >/dev/null
    makeblastdb -in sequences \
                -title "${abricate_db_name}" \
                -dbtype nucl \
                -hash_index
    popd >/dev/null
fi

################################################################################
# 5) SCAN EACH REFERENCE GENOME (USING SUMMARY TSV OR ABRicate)
################################################################################

for reference_genome in "${input_dir}"/*.{fas,fasta,fna}; do

    # Skip if no file matches
    if [[ ! -e $reference_genome ]]; then
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"

    sam_file="${output_dir}/sam/${reference_base}.sam"
    bed_file="${output_dir}/bed/${reference_base}.bed"

    # -------------------------------------------------------------------------
    # A) Use summary TSV if available and has data for this reference
    # -------------------------------------------------------------------------
    if [[ -n "${summary_tsv}" && -n "${ref2contig["${reference_name}"]}" ]]; then
        contig="${ref2contig["${reference_name}"]}"
        start0="${ref2start["${reference_name}"]}"
        endpos="${ref2end["${reference_name}"]}"

        debug_msg "Using TSV data for reference: ${reference_name}"
        debug_msg "  Contig = ${contig}, Start0 = ${start0}, End = ${endpos}"

        # Create placeholder SAM
        if [[ ! -e "${sam_file}" ]]; then
            echo -e "@HD\tVN:1.6\tSO:unsorted" > "${sam_file}"
            echo -e "@SQ\tSN:${contig}\tLN:999999999" >> "${sam_file}"
            echo -e "READ_1\t0\t${contig}\t${start0}\t60\t*\t*\t0\t0\t*\t*\tAS:i:0\tXS:i:0" \
                >> "${sam_file}"
        fi

        # Create BED file
        if [[ ! -s "${bed_file}" ]]; then
            echo -e "${contig}\t${start0}\t${endpos}\t${reference_base}\t+" \
                > "${bed_file}"
        fi

    # -------------------------------------------------------------------------
    # B) Otherwise, run ABRicate if we have no summary TSV
    # -------------------------------------------------------------------------
    else
        if [[ -z "${summary_tsv}" ]]; then
            abricate_out="${output_dir}/temp/${reference_base}_abricate.tsv"

            if [[ ! -f "${bed_file}" ]]; then
                echo "Running ABRicate on '${reference_base}'..."

                # Build the argument array for Abricate
                abricate_args=(
                    "--datadir" "${abricate_db_dir}"
                    "--db"      "${abricate_db_name}"
                    "--mincov"  "80"
                    "--minid"   "80"
                    "${reference_genome}"
                )
                
                # Run Abricate, muting stderr if debug is off
                if [[ "${debug}" -eq 1 ]]; then
                    abricate "${abricate_args[@]}" > "${abricate_out}"
                else
                    abricate "${abricate_args[@]}" > "${abricate_out}" 2>/dev/null
                fi
            fi

          # Parse ABRicate output into BED
          if [[ ! -s "${bed_file}" && -f "${abricate_out}" ]]; then
              > "${bed_file}"  # start fresh
              while IFS=$'\t' read -r line; do
                  # Skip header lines
                  [[ "$line" =~ ^# ]] && continue

                  # Example columns: FILE, SEQUENCE, START, END, GENE, ...
                  arr=($line)
                  contig="${arr[1]}"
                  startpos="${arr[2]}"
                  endpos="${arr[3]}"
                  gene="${arr[4]}"
                  strand="+"

                  echo -e "${contig}\t${startpos}\t${endpos}\t${gene}\t${strand}" \
                      >> "${bed_file}"
              done < "${abricate_out}"
          fi

          # Create placeholder SAM if needed
          if [[ ! -s "${sam_file}" ]]; then
              echo -e "@HD\tVN:1.6\tSO:unsorted" > "${sam_file}"
              if [[ -s "${bed_file}" ]]; then
                  first_line=$(head -n 1 "${bed_file}")
                  bed_contig=$(echo "$first_line" | cut -f1)
                  bed_start=$(echo "$first_line" | cut -f2)

                  echo -e "@SQ\tSN:${bed_contig}\tLN:999999999" >> "${sam_file}"
                  while IFS=$'\t' read -r c s e g str; do
                      echo -e "READ_${g}\t0\t${c}\t${s}\t60\t*\t*\t0\t0\t*\t*\tAS:i:0\tXS:i:0" \
                          >> "${sam_file}"
                  done < "${bed_file}"
              else
                  # No hits found
                  echo -e "@SQ\tSN:dummy\tLN:1" >> "${sam_file}"
              fi
          fi

        else
            debug_msg "Warning: No summary data for '${reference_name}' and no ABRicate DB provided."
            debug_msg "Skipping reference '${reference_name}'."
            continue
        fi
    fi

    # -------------------------------------------------------------------------
    # (Optional) Create genome size file for bedtools slop, etc.
    # -------------------------------------------------------------------------
    genome_size_file="${output_dir}/bed/${reference_base}_genome_size.txt"
    if [[ ! -s "${genome_size_file}" ]]; then
        awk '/^>/{if(seqname)print seqname"\t"length(seq);seqname=$1;seq="";next}
             {seq=seq$0}
             END{print seqname"\t"length(seq)}' "${reference_genome}" \
          | sed 's/>//' \
          > "${genome_size_file}"
    fi
done

################################################################################
# 6) CREATE EXTENDED BED REGIONS FOR UPSTREAM/DOWNSTREAM LD, AND GET FASTA
################################################################################

for reference_genome in "${input_dir}"/*.{fas,fasta,fna}; do
    if [[ ! -e $reference_genome ]]; then
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"

    echo "Processing extended LD for: ${reference_base}"

    bed_file="${output_dir}/bed/${reference_base}.bed"
    genome_size_file="${output_dir}/bed/${reference_base}_genome_size.txt"

    # Skip if BED file does not exist or is empty
    if [[ ! -s "${bed_file}" ]]; then
        debug_msg "  No gene found (BED empty) for '${reference_base}', skipping."
        continue
    fi

    run_num=1
    while IFS=$'\t' read -r chr start end gene_name strand; do
        query_length=$((end - start))
        debug_msg "  Creating LD for gene: ${gene_name}, #${run_num}, length: ${query_length}"

        temp_bed="${output_dir}/temp/${reference_base}_${run_num}.bed"
        echo -e "${chr}\t${start}\t${end}\t${gene_name}\t${strand}" \
            > "${temp_bed}"

        bedtools slop \
            -i "${temp_bed}" \
            -g "${genome_size_file}" \
            -b "${ld_length}" \
            >> "${output_dir}/ld_ref/${reference_base}_elongated.bed"

        run_num=$((run_num + 1))
    done < "${bed_file}"

    ld_fasta="${output_dir}/ld_ref/${reference_base}_query_gene.fas"
    bedtools getfasta \
        -fi "${reference_genome}" \
        -bed "${output_dir}/ld_ref/${reference_base}_elongated.bed" \
        -fo "${ld_fasta}"
done

################################################################################
# 7) CLEAN UP AND PREPARE LISTS FOR UNITIG-CALLER
################################################################################
echo "Removing temporary files..."
rm -rf "${output_dir}/temp"
rm -f "${input_dir}"/*.fai \
      "${input_dir}"/*.amb \
      "${input_dir}"/*.ann \
      "${input_dir}"/*.bwt \
      "${input_dir}"/*.pac \
      "${input_dir}"/*.sa

find "${output_dir}/ld_ref/" -maxdepth 1 -type f -name "*.fas" -print \
    > "${output_dir}/target_input.txt"

find "${input_dir}" -maxdepth 1 -type f \( -name "*.fas" -o -name "*.fasta" -o -name "*.fna" \) -print \
    > "${output_dir}/genome_input.txt"

################################################################################
# 8) RUN UNITIG-CALLER
################################################################################
unitig_output_pyseer="${output_dir}/unitig_output/unitig.out.pyseer"
if [[ -f "${unitig_output_pyseer}" ]]; then
    echo "unitig output '${unitig_output_pyseer}' already exists. Skipping unitig-caller."
else
    echo "Running unitig-caller..."
    unitig-caller --call \
        --refs "${output_dir}/genome_input.txt" \
        --out  "${output_dir}/unitig_output/unitig.out" \
        --kmer 31 \
        --pyseer \
        --threads "${thread_num}"
    echo "unitig-caller finished"
fi

################################################################################
# 9) CONVERT UNITIGS TO FASTA
################################################################################
unitig_fasta="${output_dir}/unitig_output/unitig.out.fasta"
awk -F'|' '{print ">" NR "\n" $1}' "${unitig_output_pyseer}" \
    > "${unitig_fasta}"

################################################################################
# 10) FILTER UNITIGS AGAINST COMBINED REFERENCE SEQUENCES USING BWA
################################################################################
bwa_output_dir="bwa_output"
mkdir -p "${output_dir}/unitig_output/${bwa_output_dir}"

combined_fasta="${output_dir}/unitig_output/combined_sequences.fasta"
cat "${output_dir}/ld_ref/"*.fas \
    > "${combined_fasta}"

bwa index "${combined_fasta}" 2>/dev/null

bwa mem -x ont2d -t "${thread_num}" \
    "${combined_fasta}" \
    "${unitig_fasta}" \
    2>/dev/null \
    > "${output_dir}/unitig_output/${bwa_output_dir}/output.sam"

bedtools bamtobed \
    -i "${output_dir}/unitig_output/${bwa_output_dir}/output.sam" \
    > "${output_dir}/unitig_output/${bwa_output_dir}/output.bed"

awk -F'\t' '!/^@/ && $3 != "*" {print $1}' \
    "${output_dir}/unitig_output/${bwa_output_dir}/output.sam" \
    | sort -n | uniq \
    > "${output_dir}/unitig_output/ids_to_extract.txt"

awk 'NR==FNR {a[$1]; next} /^>/ {header=$1; sub(/^>/, "", header); keep = !(header in a)} keep' \
    "${output_dir}/unitig_output/ids_to_extract.txt" \
    "${unitig_fasta}" \
    > "${output_dir}/unitig_output/unitig.filtered.fasta"

grep -v '^>' "${output_dir}/unitig_output/unitig.filtered.fasta" \
    > "${output_dir}/unitig_output/survived_sequences.txt"

sed -i '/^$/d' "${output_dir}/unitig_output/survived_sequences.txt"

survived_unitigs="${output_dir}/unitig_output/survived_unitigs.pyseer"
grep -F -f "${output_dir}/unitig_output/survived_sequences.txt" "${unitig_output_pyseer}" \
    > "${survived_unitigs}"

gzip "${survived_unitigs}"
rm -f "${unitig_fasta}" \
      "${output_dir}/unitig_output/unitig.filtered.fasta" \
      "${output_dir}/unitig_output/survived_sequences.txt" \
      "${output_dir}/unitig_output/ids_to_extract.txt"

echo "All tasks completed successfully."
echo "Ready to run pyseer or any further analyses."
