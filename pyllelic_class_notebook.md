---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.3
  kernelspec:
    display_name: Python 3 (ipykernel)
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

<!-- #region heading_collapsed=true -->
## Pre-setup
<!-- #endregion -->

<!-- #region heading_collapsed=true hidden=true -->
### Obtaining fastq data
<!-- #endregion -->

<!-- #region hidden=true -->
We can download rrbs (reduced representation bisulfite sequencing) data from the Encode project:
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/
<!-- #endregion -->

<!-- #region hidden=true -->
Those files are in unaligned fastq format.  We will need to align these to a reference human genome.
<!-- #endregion -->

<!-- #region hidden=true -->
### Aligning reads (using process.py)

To align reads, we'll use bowtie2 and samtools (through its pysam wrapper).

First, we need to download a genomic index sequence: http://hgdownload.soe.ucsc.edu/goldenPath/hg19
<!-- #endregion -->

```python hidden=true
# Processing imports
# from pathlib import Path
```

```python hidden=true
# Set up file paths
# index = Path(
#     "/home/andrew/allellic/hg19.p13.plusMT.no_alt_analysis_set//hg19.p13.plusMT.no_alt_analysis_set"
# )
# fastq = Path("/home/andrew/allellic/wgEncodeHaibMethylRrbsU87HaibRawDataRep1.fastq.gz")
```

<!-- #region hidden=true -->
**WARNING:** The next command is processor, RAM, and time intensive, and only needs to be run once!
<!-- #endregion -->

```python hidden=true
# Convert fastq to bam
# pyllelic.process.bowtie2_fastq_to_bam(index={bowtie_index_filename_without_suffix},
#                                       fastq={fastq_file_name},
#                                       cores=6)
```

<!-- #region hidden=true -->
Notes:

* cores is number of processor cores, adjust for your system
* instead of `out.bam` use a filename that encodes cell-line and tissue.  Our convention is: `fh_CELLLINE_TISSUE.TERT.bam`

Next, we need to sort and index the bam file using samtools functions.
<!-- #endregion -->

```python hidden=true
# Sort the bamfile
# bamfile = Path("/home/andrew/allellic/wgEncodeHaibMethylRrbsU87HaibRawDataRep1.bam")
# pyllelic.process_pysam_sort(bamfile)
```

```python hidden=true
# Create an index of the sorted bamfile
# sorted_bam = Path("")
# pyllelic.process_pysam_index(b)
```

<!-- #region hidden=true -->
Now, that sorted file (again, rename to capture cell-line and tissue info) is ready to be put in the `test` folder for analysis by pyllelic!
<!-- #endregion -->

## Set-up

```python
from pyllelic import pyllelic_class as pyllelic
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
    chrom="5",
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
data = pyllelic.GenomicPositionData(config=pyllelic.config, files_set=files_set)
```

```python
# from pathlib import PosixPath
# bam = data._bam_output[PosixPath('/home/andrew/allellic/test/fh_SW1710_URINARY_TRACT.TERT.bam')]
```

```python
# bam.values
```

```python
data.means
```

```python
data.modes
```

```python
data.diffs
```

```python
data.quma_results["SW1710"].values.head()
```

```python
# Uncomment for debugging:
", ".join(data.positions)
```

## Write Output to excel files

```python
data.save()
```

```python
# Set the filename to whatever you want
data.write_means_modes_diffs("Full1")
```

```python
import cloudpickle
```

```python
with open(r"big_data.pickle", "wb") as output_file:
    cloudpickle.dump(data, output_file)
```

```python
with open(r"big_data.pickle", "rb") as input_file:
    data2 = cloudpickle.load(input_file)
```

```python
data2.__dict__.keys()
```

```python
data2.individual_data.loc["OVK18"].dropna()["1293111"]
```

```python
large_diffs = data.diffs[(data.diffs >= 0.01).any(1)].dropna(how="all", axis=1)
```

```python
large_diffs
```

```python
temp = large_diffs.loc["HUH6"]
temp.loc[(temp != 0)].dropna()
```

```python
interesting = {}
for index, row in large_diffs.iterrows():
    if row.loc[(row != 0)].dropna().any():
        print(row.loc[(row != 0)].dropna())
        interesting[index] = row.loc[(row != 0)].dropna()
```

```python
interesting
```

```python
import pandas as pd
big_diffs = pd.DataFrame.from_dict(interesting)
```

```python
big_diffs
```

```python
big_diffs.to_excel("big_diffs.xlsx")
```

## Visualizing Data

```python
data.quma_results["SW1710"].values
```

```python
data.individual_data
```

```python
data2.histogram("OVK18", "1293111")
```

## Statistical Tests for Normality

```python
data2.summarize_allelic_data()
```

```python

```