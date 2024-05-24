import sys
import concurrent.futures
import os
from Bio import SeqIO

def generate_kmers(sequence, k):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}

def process_sequence(sequence, k):
    return generate_kmers(str(sequence.seq), k)

def process_fasta_file(file_path, k, num_threads=1):
    unique_kmers = set()

    with open(file_path, "r") as fasta_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(lambda record: process_sequence(record, k), records))

    for result in results:
        unique_kmers.update(result)

    return unique_kmers

def write_to_output_file(output_file_path, kmer_set):
    with open(output_file_path, "w") as output_file:
        for kmer in kmer_set:
            output_file.write(f"{kmer}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python your_script.py input.txt [--threads n] [--out directory/prefix]")
        sys.exit(1)

    input_file_list_path = sys.argv[1]
    kmer_size = 31  # Set the desired k-mer size

    num_threads = 1
    output_directory = os.getcwd()
    output_prefix = "output"

    # Parse command-line arguments
    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == "--threads":
            num_threads = int(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == "--out":
            output_path = sys.argv[i + 1]
            output_directory, output_prefix = os.path.split(output_path)
            i += 2
        else:
            print(f"Unrecognized option: {sys.argv[i]}")
            sys.exit(1)

    output_file_path = os.path.join(output_directory, f"{output_prefix}_output.txt")

    all_unique_kmers = set()

    with open(input_file_list_path, "r") as input_list_file:
        for fasta_file_path in input_list_file:
            fasta_file_path = fasta_file_path.strip()
            if not os.path.isfile(fasta_file_path):
                print(f"File not found: {fasta_file_path}")
                continue

            result_kmers = process_fasta_file(fasta_file_path, kmer_size, num_threads)
            all_unique_kmers.update(result_kmers)

    write_to_output_file(output_file_path, all_unique_kmers)

    print(f"All unique K-mers written to {output_file_path}")
