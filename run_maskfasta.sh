#!/bin/bash

################################################################################
# 1) PARSE COMMAND-LINE ARGUMENTS
################################################################################
while getopts q:i:d:o:t: flag
do
    case "${flag}" in
        q) query_sequence=${OPTARG};;    # Query sequence to map
        i) input_dir=${OPTARG};;         # Directory containing reference FASTA
        d) ld_length=${OPTARG};;         # Length for LD region
        o) output_dir=${OPTARG};;        # Where to store outputs
        t) thread_num=${OPTARG};;        # Number of threads
    esac
done

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
# 3) ALIGN QUERY SEQUENCE TO REFERENCES (MINIMAP2) AND CREATE BED FILES
################################################################################
for reference_genome in "${input_dir}"/*.{fas,fasta,fna}; do

    # Skip if no matching file
    if [[ ! -e $reference_genome ]]; then
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"

    # If the SAM file already exists, skip
    if [ -e "${output_dir}/sam/${reference_base}.sam" ]; then
        echo "File ${output_dir}/sam/${reference_base}.sam already exists. Skipping."
    else
        echo "Processing reference: ${reference_base}"

        # Run Minimap2 and save the SAM, suppress standard error
        minimap2 -a -t "${thread_num}" "${reference_genome}" "${query_sequence}" \
            2>/dev/null \
            > "${output_dir}/sam/${reference_base}.sam"

        # Convert SAM to BED
        bedtools bamtobed -i "${output_dir}/sam/${reference_base}.sam" \
            > "${output_dir}/bed/${reference_base}.bed"

        # Create genome length information file
        # (Chromosome name + length)
        awk '/^>/{                                \
                if (seqname)                      \
                    print seqname "\t" length(seq);\
                seqname=$1;                       \
                seq="";                           \
                next                              \
             }                                    \
             {seq = seq $0}                      \
             END {print seqname "\t" length(seq)}' "${reference_genome}" \
          | sed 's/>//; s/ / /' \
          > "${output_dir}/bed/${reference_base}_genome_size.txt"
    fi
done

################################################################################
# 4) CREATE EXTENDED BED REGIONS FOR UPSTREAM/DOWNSTREAM LD, AND GET FASTA
################################################################################
for reference_genome in "${input_dir}"/*.{fas,fasta,fna}; do

    if [[ ! -e $reference_genome ]]; then
        echo "WARNING: Fasta file for ${reference_genome} does not exist. Check the extension."
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"

    echo "Processing extended LD for: ${reference_base}"

    if [ -s "${output_dir}/bed/${reference_base}.bed" ];then

        # Counter for each gene (line in BED)
        run_num=1
    
        # Read each line (chr, start, end, gene_name, strand)
        while IFS=$'\t' read -r chr start end gene_name strand; do
            # Query gene length
            query_length=$((end - start))
            echo "  Creating LD for gene: ${gene_name}, #${run_num}, length: ${query_length}"
    
            # Temporary BED for single gene
            echo -e "${chr}\t${start}\t${end}\t${gene_name}\t${strand}" \
                > "${output_dir}/temp/${reference_base}_${run_num}.bed"
    
            # Extend regions by ld_length on each side
            bedtools slop \
                -i "${output_dir}/temp/${reference_base}_${run_num}.bed" \
                -g "${output_dir}/bed/${reference_base}_genome_size.txt" \
                -b "${ld_length}" \
                >> "${output_dir}/ld_ref/${reference_base}_elongated.bed"
    
            run_num=$((run_num + 1))
        done < "${output_dir}/bed/${reference_base}.bed"
    
        # Extract fasta for the elongated (LD) regions
        bedtools getfasta \
            -fi "${reference_genome}" \
            -bed "${output_dir}/ld_ref/${reference_base}_elongated.bed" \
            -fo "${output_dir}/ld_ref/${reference_base}_query_gene.fas"

    else
        
        echo "No gene found for ${reference_base}, Skipping" 
    fi
done

################################################################################
# 5) CLEAN UP AND PREPARE LISTS FOR UNITIG-CALLER
################################################################################
echo "Removing temporary files..."
rm -rf "${output_dir}/temp"
rm -f "${input_dir}"/*.fai \
      "${input_dir}"/*.amb \
      "${input_dir}"/*.ann \
      "${input_dir}"/*.bwt \
      "${input_dir}"/*.pac \
      "${input_dir}"/*.sa

# Find .fas files and write to target_input.txt
find "${output_dir}/ld_ref/" -maxdepth 1 -type f \( -name "*.fas" \) -print > "${output_dir}/target_input.txt"

# Find FASTA files and write to genome_input.txt
find "${input_dir}" -maxdepth 1 -type f \( -name "*.fas" -o -name "*.fasta" -o -name "*.fna" \) -print > "${output_dir}/genome_input.txt"
################################################################################
# 6) RUN UNITIG-CALLER
################################################################################

if [ -f "${output_dir}/unitig_output/unitig.out.pyseer" ];then
    echo "unitig output is exsiting, Skipping..."
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
# 7) CONVERT UNITIGS TO FASTA
################################################################################
awk -F'|' '{print ">" NR "\n" $1}' "${output_dir}/unitig_output/unitig.out.pyseer" \
    > "${output_dir}/unitig_output/unitig.out.fasta"

################################################################################
# 8) FILTER UNITIGS AGAINST COMBINED REFERENCE SEQUENCES USING BWA
################################################################################
bwa_output_dir="bwa_output"
mkdir -p "${output_dir}/unitig_output/${bwa_output_dir}"

# Combine relevant reference sequences
cat "${output_dir}/ld_ref/"*.fas \
    > "${output_dir}/unitig_output/combined_sequences.fasta"

# Index combined sequences
bwa index "${output_dir}/unitig_output/combined_sequences.fasta" 2>/dev/null

# Align unitigs to combined sequences, suppress BWA printing
bwa mem -k 13 -t 6 \
    "${output_dir}/unitig_output/combined_sequences.fasta" \
    "${output_dir}/unitig_output/unitig.out.fasta" \
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
    "${output_dir}/unitig_output/unitig.out.fasta" \
    > "${output_dir}/unitig_output/unitig.filtered.fasta"

# Extract only the sequences (not headers)
grep -v '^>' "${output_dir}/unitig_output/unitig.filtered.fasta" \
    > "${output_dir}/unitig_output/survived_sequences.txt"

# Remove empty lines
sed -i '/^$/d' "${output_dir}/unitig_output/survived_sequences.txt"

################################################################################
# 9) GREP CORRESPONDING LINES FROM PYSEER UNITIG OUTPUT
################################################################################
grep -F -f "${output_dir}/unitig_output/survived_sequences.txt" \
    "${output_dir}/unitig_output/unitig.out.pyseer" \
    > "${output_dir}/unitig_output/survived_unitigs.pyseer"

# Gzip the final set
if [ -e "${output_dir}/unitig_output/survived_unitigs.pyseer.gz" ]; then
    rm "${output_dir}/unitig_output/survived_unitigs.pyseer.gz"
fi
gzip "${output_dir}/unitig_output/survived_unitigs.pyseer"

################################################################################
# 10) FINAL CLEANUP
################################################################################
echo "Remove temporary unitig files..."
rm -f "${output_dir}/unitig_output/unitig.out.fasta" \
      "${output_dir}/unitig_output/unitig.filtered.fasta" \
      "${output_dir}/unitig_output/survived_sequences.txt" \
      "${output_dir}/unitig_output/ids_to_extract.txt"

echo "Ready to run pyseer or any further analyses."
