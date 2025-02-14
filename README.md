# maskGWAS

**maskGWAS** is a add-on tool of [unitig-caller](https://github.com/bacpop/unitig-caller) and [pyseer](https://github.com/weecology/pyseer) designed for Genome-Wide Epistatic Studies (GWES) by masking target genes in GWAS. It generates unitig files of whole genomes by masking the specified target gene. The primary goal is to detect hidden epistatic relationships involving the target gene.

## Warning

1. **Current Compatibility**: `maskGWAS` is currently only runnable on if the requirements are fully ready. Check the requirements first, if there is any error.
2. **Dependency Requirement**: Ensure that `unitig-caller` is properly downloaded and configured in your environment.

## Requirements

Before using `maskGWAS`, ensure the following software is installed:

1. **[minimap2](https://github.com/lh3/minimap2)** (version 2.24 or higher)
2. **[BEDtools](https://bedtools.readthedocs.io/en/latest/)** (version 2.30.0 or higher)
3. **[PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)** (Note: PLINK 2 is **not** supported)
4. **[unitig-caller](https://github.com/bacpop/unitig-caller)** (or higher)
5. **[BWA](http://bio-bwa.sourceforge.net/)** (version 0.7.17 or higher)
6. **[pyseer](https://github.com/weecology/pyseer)** (version 1.3.11 or higher)
7. **[ABRICATE](https://github.com/tseemann/abricate)** (version 1.0.0 or higher)

## Installation

Clone the repository to your local directory:

```ruby
git clone https://github.com/jeju2486/maskGWAS.git
cd maskGWAS
```

## Usage

### 1. Calculate Linkage Disequilibrium (LD)

#### a. Using Fasta Files

To calculate LD using fasta files in your directory, run:

```ruby
bash run_ld_calculation.sh -i "$input_dir" -o "$output_dir" -t $SLURM_CPUS_PER_TASK
```

#### b. Using Core Alignment Files (Optional)

If you have a core alignment file from PIRATE or any other source, you can reduce variant calling time:

```ruby
bash run_ld_calculation_from_aln.sh -i "$pirate_result_dir" -o "$output_dir" -t $SLURM_CPUS_PER_TASK
```

### 2.Visualize LD Decay (Optional)

Visualize LD decay using `ggplot2` in R:

```ruby
Rscript plotting_lddecay.r -i "$output_dir/ld_results_sampled.ld" -o "$output_dir"
```
**Notes:**
- This generates LD decay images.
- Ensure the first file in the directory is your reference genome.
- *Todo*:
  - Make the reference genome selectable.
  - Add debugging messages to verify environment setup.


### 3.Mask Target Gene
Run the masking script to generate SAM/BED files and LD information:

```ruby
bash run_maskfasta.sh \
  -g "/path/to/query_sequence.fasta" \
  -i "/path/to/input_dir" \
  -d 3000 \
  -o "/path/to/output_dir" \
  -t $SLURM_CPUS_PER_TASK
```

#### If there is specific region you want to mask (Optional) 

If you have the specific region you want to mask, you can use the `-s` parameter and give the file sepecify the region. The file should be tab-delimited format and should specify the contig, starting and end point. Like:

`target region.tsv`

```
Reference	Gene_presence	Contig	Start0	End
GCF_000144955.fas	Yes	NC_017338.2	33707	57871
GCF_000159535.fas	Yes	NC_017342.1	459032	486212
GCF_000189455.fas	No	NZ_CP025395.1	33931	34407
GCF_000568455.fas	Yes	NZ_CP007176.1	34275	88793
```

Then the running code will be like:

```ruby
bash run_maskfasta.sh\
  -s "path/to/target_region.tsv" \
  -i "/path/to/input_dir" \
  -d 3000 \
  -o "/path/to/output_dir" \
  -t $SLURM_CPUS_PER_TASK
```

**Output Directories:**

- **`./sam/`**  
  Contains Sequence Alignment/Map (SAM) files.
- **`./bed/`**  
  Contains BED files (simplified SAM files).
- **`./ld_ref/`**  
  Contains LD information files, including BED and FASTA files for LD regions.
- **`./unitig_output/`**  
  Contains unitig outputs:
  - `unitig.output.pyseer`: Unitigs of whole genomes for analysis.
  - `survived_unitigs.pyseer.gz`: Unitig files without LD regions, converted for pyseer.

### 4. Run pyseer

Execute pyseer and sort the result files, separating extreme outliers for plotting:

```ruby
bash run_pyseer.sh \
  -t "$treefile" \
  -p "$phenotype" \
  -i "$maskfasta_output_dir/unitig_output" \
  -o "$pyseer_output_dir" \
  -P "prefix" \          # Optional
  -T "$SLURM_CPUS_PER_TASK" \  # Optional
  -s "$pyseer_script"
```

#### Alternative: Use `run_code.sh`

You can use the all-in-one `run_code.sh` script, which includes instructions for required program versions.  
**Warning**: Remember to update file names accordingly.

```ruby
bash run_code.sh
```

### 5. Plot pyseer Results

**Note**: This step is currently only runnable on ARC servers. Ensure you update directory paths accordingly.

```ruby
bash run_pyseer_plotting.sh \
  -r "$reference_dir" \
  -g "$gff_dir" \
  -s "$pyseer_script" \
  -o "$pyseer_output_dir" \
  -p "sccmec"  # Optional: Prefix for output files
```

**Notes:**
- The prefix should ideally match the one used in `run_pyseer.sh`.
- Generates two plots: one for extreme outliers and one for normal data.
- *Todo*: Make the plotting process more adaptable.

## Output

- **`./sam/`**  
  Sequence Alignment/Map (SAM) files.
- **`./bed/`**  
  BED files (simplified SAM files).
- **`./ld_ref/`**  
  LD information files, including BED and FASTA files for LD regions.
- **`./unitig_output/`**  
  1. `unitig.output.pyseer`: Unitigs of whole genomes used for analysis.
  2. `potential.kmer_output.txt`: k-mers from LD regions.
  3. `extracted.unitig.out.pyseer`: Unitig files without LD regions, formatted for pyseer.

## Running PYSEER

`maskGWAS` automates running pyseer and sorting results, separating extreme outliers for individual plots.

```ruby
bash run_pyseer.sh \
  -t "$treefile" \
  -p "$phenotype" \
  -i "$maskfasta_output_dir/unitig_output" \
  -o "$pyseer_output_dir" \
  -P "prefix" \          # Optional
  -T "$SLURM_CPUS_PER_TASK" \  # Optional
  -s "$pyseer_script"
```

## Plotting PYSEER Results

Generate plots for pyseer results. Ensure directory paths are correctly set.

```ruby
bash run_pyseer_plotting.sh \
  -r "$input_fasta_dir" \
  -g "$gff_dir" \
  -s "$pyseer_script" \
  -o "$pyseer_output_dir" \
  -p "sccmec"  # Optional: Prefix for output files
```
This code will automatically find the gff files in the given directory corresponding to the fasta file and assuming all of them are reference genome but not the sample one (Read the [pyseer tutorial page](https://pyseer.readthedocs.io/en/master/tutorial.html#k-mer-association-with-mixed-effects-model) for more detail)

**Output:**
- Two plots:
  1. Extreme outliers.
  2. Normal distribution.

*Todo*: Enhance plot customization and adaptability.

## To-Do List

- **Generalization**: Learn Conda to make `maskGWAS` runnable outside ARC servers.
- **Reference Genome Selection**: Allow users to select the reference genome file.
- **Plotting Enhancements**: Make plotting scripts more adaptable and customizable.

## Additional Information

- **Repository**: [maskGWAS on GitHub](https://github.com/jeju2486/maskGWAS)
- **Main Tool**: Utilizes [unitig-caller](https://github.com/bacpop/unitig-caller) as the primary tool for generating unitig files.

## Contact

For questions, issues, or contributions, please open an issue on the [GitHub repository](https://github.com/jeju2486/maskGWAS) or contact the maintainer.
