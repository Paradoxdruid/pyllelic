---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Example Pyllelic Use-Case Notebook

<!-- #region heading_collapsed=true -->
## Background
<!-- #endregion -->

<!-- #region hidden=true -->
This notebook illustrates the import and use of `pyllelic` in a jupyter environment.

See https://github.com/Paradoxdruid/pyllelic for further details.
<!-- #endregion -->

## Pre-setup

<!-- #region heading_collapsed=true -->
### Obtaining fastq data
<!-- #endregion -->

<!-- #region hidden=true -->
We can download rrbs (reduced representation bisulfite sequencing) data from the Encode project:
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/
<!-- #endregion -->

<!-- #region hidden=true -->
Those files are in unaligned fastq format.  We will need to align these to a reference human genome.
<!-- #endregion -->

### Aligning reads (using process.py)

To align reads, we'll use bowtie2 and samtools (through its pysam wrapper).

First, we need to download a genomic index sequence: http://hgdownload.soe.ucsc.edu/goldenPath/hg19

```python
# Processing imports
# from pathlib import Path
```

```python
# Set up file paths
# index = Path(
#     "/home/andrew/allellic/hg19.p13.plusMT.no_alt_analysis_set//hg19.p13.plusMT.no_alt_analysis_set"
# )
# fastq = Path("/home/andrew/allellic/wgEncodeHaibMethylRrbsU87HaibRawDataRep1.fastq.gz")
```

**WARNING:** The next command is processor, RAM, and time intensive, and only needs to be run once!

```python
# Convert fastq to bam
# pyllelic.process.bowtie2_fastq_to_bam(index={bowtie_index_filename_without_suffix},
#                                       fastq={fastq_file_name},
#                                       cores=6)
```

Notes:

* cores is number of processor cores, adjust for your system
* instead of `out.bam` use a filename that encodes cell-line and tissue.  Our convention is: `fh_CELLLINE_TISSUE.TERT.bam`

Next, we need to sort and index the bam file using samtools functions.

```python
# Sort the bamfile
# bamfile = Path("/home/andrew/allellic/wgEncodeHaibMethylRrbsU87HaibRawDataRep1.bam")
# pyllelic.process_pysam_sort(bamfile)
```

```python
# Create an index of the sorted bamfile
# sorted_bam = Path("")
# pyllelic.process_pysam_index(b)
```

Now, that sorted file (again, rename to capture cell-line and tissue info) is ready to be put in the `test` folder for analysis by pyllelic!

## Set-up

```python
import pyllelic
```

```python
# set up your disk location:
# base_path should be the directory we'll do our work in
# make a sub-directory under base_path with a folder named "test"
# and put the .bam and .bai files in "test"

# OSX setup
# pyllelic.set_up_env_variables(
#     base_path="/Users/abonham/documents/test_allelic/",
#     prom_file="TERT-promoter-genomic-sequence.txt",
#     prom_start="1293000",
#     prom_end="1296000",
#     chrom="5",
#     offset=1298163,
# )

# Windows set-up
pyllelic.set_up_env_variables(
    base_path="/home/andrew/allellic/",
    prom_file="TERT-promoter-genomic-sequence.txt",
    prom_start="1293000",
    prom_end="1296000",
    chrom="chr5",
    offset=1298163,
)
```

## Main Parsing Functions

```python
files_set = pyllelic.make_list_of_bam_files()  # finds bam files
```

```python
# Uncomment for debugging:
files_set
```

```python
# index bam and creates bam_output folders/files
# set process False to skip writing output files if they already exist
positions = pyllelic.index_and_fetch(files_set, process=True)
```

```python
# Turn off pretty printing to see position list better
%pprint
```

```python
# Uncomment for debugging:
positions
```

```python
# Turn back on pretty printing
%pprint
```

```python
# Only needs to be run once, generates static files
# pyllelic.genome_parsing()

# Can also take sub-list of directories to process
# pyllelic.genome_parsing([pyllelic.config.bam_directory / "fh_BONHAM_TISSUE.TERT.bam"])
```

```python
cell_types = pyllelic.extract_cell_types(files_set)
```

```python
# Uncomment for debugging
cell_types
```

```python
# Set filename to whatever you want
df_list = pyllelic.run_quma_and_compile_list_of_df(
    cell_types, "test2.xlsx",
    run_quma=True,
)  # to skip quma: , run_quma=False)
```

```python
# Uncomment for debugging
df_list.keys()
```

```python
df_list["SORTED"]
```

```python
means = pyllelic.process_means(df_list, positions, files_set)
```

```python
# Uncomment for debugging
means
```

```python
modes = pyllelic.process_modes(df_list, positions, files_set)
```

```python
# Uncomment for debugging
modes
```

```python
diff = pyllelic.find_diffs(means, modes)
```

```python
# Uncomment for debugging
diff
```

## Write Output to excel files

```python
# Set the filename to whatever you want
pyllelic.write_means_modes_diffs(means, modes, diff, "Test1")
```

## Visualizing Data

```python
final_data = pyllelic.pd.read_excel(
    pyllelic.config.base_directory.joinpath("Test1_diff.xlsx"),
    dtype=str,
    index_col=0,
    engine="openpyxl",
)
```

```python
final_data
```

```python
individual_data = pyllelic.return_individual_data(df_list, positions, files_set)
```

```python
# Uncomment for debugging
individual_data
```

```python
individual_data.loc["NCIH196", "1295937"]
```

```python
individual_data.loc["SORTED"]["1295680"]
```

```python
individual_data.loc["CALU1"]
```

```python
pyllelic.histogram(individual_data, "SORTED", "1295680")
```

```python
pyllelic.histogram(individual_data, "SORTED", "1295903")
```

```python
pyllelic.histogram(individual_data, "SW1710", "1295089")
```

```python
pyllelic.histogram(individual_data, "CALU1", "1295937")
```

```python
pyllelic.histogram(individual_data, "NCIH196", "1295937")
```

```python
pyllelic.histogram(individual_data, "NCIH196", "1294945")
```

```python
final_data.loc["SW1710"]
```

## Statistical Tests for Normality

```python
pyllelic.summarize_allelic_data(individual_data, diff)
```

```python

```
