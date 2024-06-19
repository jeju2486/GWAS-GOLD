# maskGWAS
GWES by making target gene.

This will generate the unitig file of a whole genome, which can be used for pyseer, with masking the target gene. This uses [unitig-caller](https://github.com/bacpop/unitig-caller) as a main tool. Main goal is detecting the hidden epistatic relationship of target gene.

warning
1. This is currently only runable under ARC server (#todo-list learn conda to make it runable generally) 
2. This assumes there is unitig-caller downloaded properly under your environment. 

# Usage

First, you need to clone the whole respository into your local directory

```ruby
git clone https://github.com/jeju2486/maskGWAS.git
```

Then run the `run_plink` file to calculate LD (linkage disequilibrium) length. This will detect fasta files in your directory.

```ruby
bash run_plink.sh -i "/path/to/input_dir" -o "/path/to/output_dir"
```

This will generate the image files of ld decaying. You need to choose the proper threshold.This assumes the first file in the directory as reference genome to run it. (#todo-list make it changable). This assumes you loaded the unitig-caller and conda environemt properly before running it (#todo-list: add the debugging message to check if they really loaded properly).

```ruby
bash run_maskfasta.sh -q "/path/to/query_sequence.fasta" -i "/path/to/input_dir" -d 3000 -o "/path/to/output_dir" -t 12
```

This will generate two directores for sam/bed files and one directory (ld_ref) for ld information files

This also generate unitig_output file whch containsg the unitig informations. You need unitig_output/extracted.unitig.out.pyseer file for pyseer running. 

# Output

* `**./sam**`

folder of Sequce Alignment/Map (SAM) files 

* `**./bed/**`

folder of BED files. BED file is the simplified file of SAM file

* `**./ld_ref/**`

folder containing ld information. It includes the bed and fasata files of LD (linkage disequilibrium) sequences

* `**./unitig_output/**`

folder containing unitig outputs. It includes

1. `unitig.output.pyseer`

unitig of whole genomes used for anaylsis

2. `potential.kmer_output.txt`

k-mers of LD region

3. `extracted.unitig.out.pyseer`

unitig file of genome without LD region. converted into pyseer format.

# For running PYSEER

The code below are for helping you to run the pyseer. Use for reference only. For more detail read following [link](https://pyseer.readthedocs.io/en/master/tutorial.html).


```ruby
pyseer --lmm --phenotypes "$phenotype" --kmers "$unitig_result"/extracted.unitig.out.pyseer.gz --similarity phylogeny_wholegenome.tsv --print-samples --output-patterns "$result_dir"/kmer_patterns.txt --cpu 24 > "$result_dir"/sccmec_kmers.txt
```
