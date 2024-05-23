# maskGWAS
GWES by making target gene

warning: it is currently only runable under ARC server (#todo-list learn conda to make it runable generally)

# usage

```ruby
bash run_plink.sh -i <isolate_dir> -o <output_dir>
```

This will generate the image files of ld decaying. You need to choose the proper threshold.This assumes the first file in the directory as reference file to run it. (#todo-list make it changable)

```ruby
bash run_maskfasta.sh -q "/path/to/query_sequence.fasta" -i "/path/to/input_dir" -d 3000 -o "/path/to/output_dir"
```


