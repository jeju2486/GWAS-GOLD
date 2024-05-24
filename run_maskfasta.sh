#!/bin/bash

# Load necessary modules
module purge
module load minimap2/2.24-GCCcore-11.3.0
module load BEDTools/2.30.0-GCC-11.2.0
module load Anaconda3
source activate $DATA/myenv

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

# Create necessary directories
mkdir -p "$output_dir/sam" "$output_dir/bed" "$output_dir/ld_ref" "$output_dir/temp"

rm "$output_dir/ld_ref/*"

# Iterate through reference genomes
for reference_genome in "$input_dir"/*.{fas,fasta,fna}; do
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
        minimap2 -a "$reference_genome" "$query_sequence" > "$output_dir/sam/${reference_base}.sam"
        
        # Convert SAM to BED
        bedtools bamtobed -i "$output_dir/sam/${reference_base}.sam" > "$output_dir/bed/${reference_base}.bed"

        # Create genome length information file
        awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' "$reference_genome" | sed 's/>//; s/ / /' > "$output_dir/bed/${reference_base}_genome_size.txt"
    fi
done

# Create BED files for upstream_ld, downstream_ld, and query gene itself
for reference_genome in "$input_dir"/*.{fas,fasta,fna}; do 
    # Check if the file exists to handle the case when there are no matching files
    if [[ ! -e $reference_genome ]]; then
        echo "Fasta file for ${reference_genome} does not exist. Check the extension. Ignore if this is not your extension."
        continue
    fi

    reference_name=$(basename "$reference_genome")
    reference_base="${reference_name%.*}"
    
    echo "Processing ${reference_base} ..."
    # Count the run
    run_num=1
    
    # Iterate over each line in the BED file
    while IFS=$'\t' read -r chr start end gene_name strand; do

        # Calculate the query gene length
        query_length=$((end - start))
        
        echo "Processing gene: ${gene_name}, run = ${run_num}, length: ${query_length}"

        # Create temporary BED file for the current gene
        echo -e "${chr}\t${start}\t${end}\t${gene_name}\t${strand}" > "$output_dir/temp/${reference_base}_${run_num}.bed"

        # Create upstream_ld and downstream_ld BED files, adjusting for the query gene length
        bedtools slop -i "$output_dir/temp/${reference_base}_${run_num}.bed" -g "$output_dir/bed/${reference_base}_genome_size.txt" -l "${ld_length}" -r $((-query_length)) -s >> "$output_dir/ld_ref/${reference_base}_upstream_ld.bed"
        bedtools slop -i "$output_dir/temp/${reference_base}_${run_num}.bed" -g "$output_dir/bed/${reference_base}_genome_size.txt" -l $((-query_length)) -r "${ld_length}" -s >> "$output_dir/ld_ref/${reference_base}_downstream_ld.bed"
        
        # Increment the run number
        run_num=$((run_num + 1))
        
    done < "$output_dir/bed/${reference_base}.bed"
    
    # Get FASTA sequences for the query gene and ld regions
    bedtools getfasta -fi "$reference_genome" -bed "$output_dir/bed/${reference_base}.bed" -fo "$output_dir/ld_ref/${reference_base}_query_gene.fas"
    bedtools getfasta -fi "$reference_genome" -bed "$output_dir/ld_ref/${reference_base}_upstream_ld.bed" -fo "$output_dir/ld_ref/${reference_base}_upstream_ld.fas"
    bedtools getfasta -fi "$reference_genome" -bed "$output_dir/ld_ref/${reference_base}_downstream_ld.bed" -fo "$output_dir/ld_ref/${reference_base}_downstream_ld.fas"
        
done

# Clean up temporary files
rm -r "$output_dir/temp"
rm "$output_dir"/*txt

echo "Masking Finished, Now run unitig-caller..."

#get fasta file addresses from input directory
ls -d -1 "$output_dir"/ld_ref/*.fas > "$output_dir"/target_input.txt
ls -d -1 "$input_dir"/*.{fas,fasta,fna} > "$output_dir"/genome_input.txt

mkdir -p "$output_dir/unitig_output"

echo "Running unitig-caller for genomes"

#run unitig-caller
unitig-caller --call --refs "$output_dir"/genome_input.txt --out "$output_dir/unitig_output"/unitig.out --kmer 31 --pyseer --threads "${thread_num}"

#run python code
python generate_kmer.py "$output_dir"/target_input.txt --threads "${thread_num}" --out "$output_dir/unitig_output"/potential.kmers

#extract the common field between the refernece genome and sccmec k-mers and delete the row from reference unitig result regardless of the second field
grep -vf "$output_dir/unitig_output"/potential.kmers_output.txt "$output_dir/unitig_output"/unitig.out.pyseer > "$output_dir/unitig_output"/extracted.unitig.out.pyseer

#gzip the file
gzip "$output_dir/unitig_output"/extracted.unitig.out.pyseer

echo "unitig-caller finished"

echo "Ready to run pyseer"
