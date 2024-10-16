#!/bin/bash

# Enable nullglob to handle cases where no files match
shopt -s nullglob

# Load necessary modules
module purge
module load minimap2/2.24-GCCcore-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0
module load Anaconda3/2023.09-0
source activate /data/biol-micro-genomics/kell7366/python3_12_2

# Parse command-line arguments
while getopts q:i:d:o:t: flag
do
    case "${flag}" in
        q) query_sequence=${OPTARG};;
        i) input_dir=${OPTARG};;
        d) ld_length=${OPTARG};;
        o) output_dir=${OPTARG};;
        t) thread_num=${OPTARG};;
    esac
done

# Validate command-line arguments
if [[ -z "$query_sequence" || -z "$input_dir" || -z "$ld_length" || -z "$output_dir" || -z "$thread_num" ]]; then
    echo "Usage: $0 -q query_sequence -i input_dir -d ld_length -o output_dir -t thread_num"
    exit 1
fi

if ! [[ "$ld_length" =~ ^[0-9]+$ ]]; then
    echo "Error: ld_length must be a positive integer."
    exit 1
fi

# Create necessary directories
mkdir -p "$output_dir/sam" "$output_dir/bed" "$output_dir/ld_ref" "$output_dir/temp"

# Remove existing files in ld_ref directory
rm -f "$output_dir/ld_ref"/*

# Iterate through reference genomes
for extension in fas fasta fna; do
    for reference_genome in "$input_dir"/*.$extension; do
        # Check if the file exists to handle the case when there are no matching files
        if [[ ! -e $reference_genome ]]; then
            continue
        fi

        reference_name=$(basename "$reference_genome")
        reference_base="${reference_name%.*}"
        
        if [ -e "$output_dir/sam/${reference_base}.sam" ]; then
            echo "File $output_dir/sam/${reference_base}.sam already exists. Skipping."
        else
            echo "Processing $reference_base"

            # Run Minimap2 and save the SAM file
            minimap2 -a "$reference_genome" "$query_sequence" -t "$thread_num" > "$output_dir/sam/${reference_base}.sam"
            
            # Convert SAM to BED
            bedtools bamtobed -i "$output_dir/sam/${reference_base}.sam" > "$output_dir/bed/${reference_base}.bed"

            # Create genome length information file
            awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' "$reference_genome" | sed 's/>//; s/ / /' > "$output_dir/bed/${reference_base}_genome_size.txt"
        fi
    done
done

# Create BED files for upstream_ld, downstream_ld, and query gene itself
for extension in fas fasta fna; do
    for reference_genome in "$input_dir"/*.$extension; do
        # Check if the file exists to handle the case when there are no matching files
        if [[ ! -e $reference_genome ]]; then
            echo "Fasta file for ${reference_genome} does not exist. Check the extension. Ignore if this is not your extension."
            continue
        fi

        reference_name=$(basename "$reference_genome")
        reference_base="${reference_name%.*}"
        
        echo "Processing ${reference_base} ..."
        # Initialize run number
        run_num=1
        
        # Iterate over each line in the BED file
        while IFS=$'\t' read -r chr start end; do
            gene_name="${reference_base}"
            strand="."
            # Calculate the query gene length
            query_length=$((end - start))
            
            echo "Processing gene: ${gene_name}, run = ${run_num}, length: ${query_length}"

            # Create temporary BED file for the current gene
            echo -e "${chr}\t${start}\t${end}\t${gene_name}\t${strand}" > "$output_dir/temp/${reference_base}_${run_num}.bed"

            # Create upstream_ld and downstream_ld BED files, adjusting for the query gene length
            bedtools slop -i "$output_dir/temp/${reference_base}_${run_num}.bed" -g "$output_dir/bed/${reference_base}_genome_size.txt" -b "${ld_length}" >> "$output_dir/ld_ref/${reference_base}_elongated.bed"
            
            # Increment the run number
            run_num=$((run_num + 1))
            
        done < "$output_dir/bed/${reference_base}.bed"
        
        # Get FASTA sequences for the query gene and ld regions
        bedtools getfasta -fi "$reference_genome" -bed "$output_dir/ld_ref/${reference_base}_elongated.bed" -fo "$output_dir/ld_ref/${reference_base}_query_gene.fas"       
    done
done

echo "Removing temporary files"

# Clean up temporary files
rm -r "$output_dir/temp"
rm -f "$input_dir"/*.fai "$input_dir"/*.amb "$input_dir"/*.ann "$input_dir"/*.bwt "$input_dir"/*.pac "$input_dir"/*.sa

# Activate the next Anaconda environment
module purge
module load Anaconda3/2024.02-1
source activate /data/biol-micro-genomics/kell7366/python3_12_2

# Get fasta file addresses from input directories
ls -d -1 "$output_dir"/ld_ref/*.fas > "$output_dir"/target_input.txt
ls -d -1 "$input_dir"/*.{fas,fasta,fna,fa} > "$output_dir"/genome_input.txt

# Create unitig_output directory
mkdir -p "$output_dir/unitig_output"

echo "Running unitig-caller for genomes"

# Run unitig-caller
unitig-caller --call --refs "$output_dir/genome_input.txt" --out "$output_dir/unitig_output/unitig.out" --kmer 31 --pyseer --threads "$thread_num"

echo "unitig-caller finished"
echo "Processing unitigs"

# Load modules for unitig processing
module purge
module load BWA/0.7.17-GCCcore-11.2.0
module load BEDTools/2.30.0-GCC-11.2.0
module load SAMtools/1.16.1-GCC-11.3.0

# Convert unitigs to fasta format
awk -F'|' '{print ">" NR "\n" $1}' "$output_dir/unitig_output/unitig.out.pyseer" > "$output_dir/unitig_output/unitig.out.fasta"

echo "Filtering out the unitigs"

# Create an output directory for BWA results
bwa_output_dir="$output_dir/unitig_output/bwa_output"
mkdir -p "$bwa_output_dir"

# Combine sequences from the ld_ref directory
cat "$output_dir/ld_ref/"*.fas > "$output_dir/combined_sequences.fasta"

# Index the combined sequences
bwa index "$output_dir/combined_sequences.fasta"

# Run bwa mem
bwa mem -k 20 -t 6 "$output_dir/combined_sequences.fasta" "$output_dir/unitig_output/unitig.out.fasta" > "${bwa_output_dir}/output.sam"

# Filter SAM file for exact matches
samtools view -h "${bwa_output_dir}/output.sam" | awk '/^@/ || $6 ~ /^[0-9]+M$/' > "${bwa_output_dir}/exact_matches.sam"

# Convert filtered SAM to BED
bedtools bamtobed -i "${bwa_output_dir}/exact_matches.sam" > "${bwa_output_dir}/output.bed"

echo "Extracting ids"

# Extract unique IDs from exact_matches.sam
awk -F'\t' '!/^@/ && $3 != "*" {print $1}' "${bwa_output_dir}/exact_matches.sam" | sort -n | uniq > "$output_dir/unitig_output/ids_to_extract.txt"

# Filter unitig.out.fasta to exclude IDs in ids_to_extract.txt
awk 'NR==FNR {a[$1]; next} /^>/ {header=$1; sub(/^>/, "", header); keep = !(header in a)} keep' "$output_dir/unitig_output/ids_to_extract.txt" "$output_dir/unitig_output/unitig.out.fasta" > "$output_dir/unitig_output/unitig.filtered.fasta"

# Extract sequences from unitig.filtered.fasta
grep -v '^>' "$output_dir/unitig_output/unitig.filtered.fasta" > "$output_dir/unitig_output/survived_sequences.txt"
sed -i '/^$/d' "$output_dir/unitig_output/survived_sequences.txt"

# Extract corresponding lines from unitig.out.pyseer
grep -F -f "$output_dir/unitig_output/survived_sequences.txt" "$output_dir/unitig_output/unitig.out.pyseer" > "$output_dir/unitig_output/survived_unitigs.pyseer"

# Gzip the survived_unitigs.pyseer file
if [ -e "$output_dir/unitig_output/survived_unitigs.pyseer.gz" ]; then
    rm "$output_dir/unitig_output/survived_unitigs.pyseer.gz"
fi
gzip "$output_dir/unitig_output/survived_unitigs.pyseer"

# Replace the original unitig.out.pyseer.gz with the filtered version
mv "$output_dir/unitig_output/survived_unitigs.pyseer.gz" "$output_dir/unitig_output/unitig.out.pyseer.gz"

echo "Removing temporary files"
rm "$output_dir/unitig_output/unitig.out.fasta" "$output_dir/unitig_output/unitig.filtered.fasta" "$output_dir/unitig_output/survived_sequences.txt" "$output_dir/unitig_output/ids_to_extract.txt"
rm "$output_dir/combined_sequences.fasta"
rm -r "$bwa_output_dir"

echo "Unitig filtering completed"

# Proceed to run pyseer
echo "Ready to run pyseer"
