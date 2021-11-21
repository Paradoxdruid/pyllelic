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

Source code: <https://github.com/Paradoxdruid/pyllelic>

Documentation: <https://paradoxdruid.github.io/pyllelic/>
<!-- #endregion -->

<!-- #region heading_collapsed=true -->
## Pre-setup / File preparation
<!-- #endregion -->

<!-- #region heading_collapsed=true hidden=true -->
### Obtaining fastq data
<!-- #endregion -->

<!-- #region hidden=true -->
We can download rrbs (reduced representation bisulfite sequencing) data from the Encode project:
<http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/>
<!-- #endregion -->

<!-- #region hidden=true -->
Those files are in unaligned fastq format.  We will need to align these to a reference human genome.
<!-- #endregion -->

<!-- #region heading_collapsed=true hidden=true -->
### Aligning reads (using process.py)

To align reads, we'll use bowtie2 and samtools (through its pysam wrapper).

First, we need to download a genomic index sequence: <http://hgdownload.soe.ucsc.edu/goldenPath/hg19>
<!-- #endregion -->

<!-- #region hidden=true -->
```python
# Processing imports
from pathlib import Path
```
<!-- #endregion -->

<!-- #region hidden=true -->
```python
# Set up file paths
index = Path(
    "/home/andrew/allellic/hg19.p13.plusMT.no_alt_analysis_set//hg19.p13.plusMT.no_alt_analysis_set"
)
fastq = Path("/home/andrew/allellic/wgEncodeHaibMethylRrbsU87HaibRawDataRep1.fastq.gz")
```
<!-- #endregion -->

<!-- #region hidden=true -->
**WARNING:** The next command is processor, RAM, and time intensive, and only needs to be run once!
<!-- #endregion -->

<!-- #region hidden=true -->
```python
# Convert fastq to bam
from pyllelic import process
process.bowtie2_fastq_to_bam(index={bowtie_index_filename_without_suffix},
                                      fastq={fastq_file_name},
                                      cores=6)
```
<!-- #endregion -->

<!-- #region hidden=true -->
Notes:

* cores is number of processor cores, adjust for your system
* instead of `out.bam` use a filename that encodes cell-line and tissue.  Our convention is: `fh_CELLLINE_TISSUE.TERT.bam`

Next, we need to sort and index the bam file using samtools functions.
<!-- #endregion -->

<!-- #region hidden=true -->
```python
# Sort the bamfile
bamfile = Path("/home/andrew/allellic/wgEncodeHaibMethylRrbsU87HaibRawDataRep1.bam")
process.process_pysam_sort(bamfile)
```
<!-- #endregion -->

<!-- #region hidden=true -->
```python
# Create an index of the sorted bamfile
sorted_bam = Path("")
process.process_pysam_index(b)
```
<!-- #endregion -->

<!-- #region hidden=true -->
Now, that sorted file (again, rename to capture cell-line and tissue info) is ready to be put in the `test` folder for analysis by pyllelic!
<!-- #endregion -->

## Set-up

```python
from pyllelic import pyllelic
```

1. Set up your disk location:  ```base_path``` should be the directory we'll do our work in
2. Make a sub-directory under ```base_path``` with a folder named ```test``` and put the ```.bam``` and ```.bai``` files in ```test```

```python
config = pyllelic.set_up_env_variables(
    base_path="/home/andrew/allellic/",  # Windows WSL set-up
#     base_path="/Users/abonham/documents/test_allelic/",  # OSX setup
    prom_file="TERT-promoter-genomic-sequence.txt",
    prom_start="1293000",
    prom_end="1296000",
    chrom="5",
    offset=1298163,
)
```

## Main Parsing Functions

### Find files to analyze

```python
files_set = pyllelic.make_list_of_bam_files(config)
```

```python
# files_set # uncomment for debugging
```

```python
files_set = files_set[0:2]  # grab only first two files for quick run
```

### Perform full methylation analysis and generate data object

```python
data = pyllelic.GenomicPositionData(config=config, files_set=files_set)
```

### Check main data outcomes

```python
data.means.head()
```

```python
", ".join(data.positions)
```

```python
data.modes.head()
```

```python
data.diffs.head()
```

```python
data.quma_results[data.means.index[0]].values.head()
```

## Write Output

### Save entire object as pickle

```python
# data.save_pickle("test3_data.pickle")
```

### Save Quma Results to excel

```python
# data.save("output.xlsx")
```

### Save analysis files (means, modes, diffs) to excel

```python
# data.write_means_modes_diffs("Full1") # Sets the filename stub
```

### Reopen saved object

```python
# data = pyllelic.GenomicPositionData.from_pickle("test3_data.pickle")
```

```python
# old import system, do not use
# with open(r"big_data.pickle", "rb") as input_file:
#     data = cloudpickle.load(input_file)
```

## Initial Data Analysis

### View raw data of methylation percentage per read

```python
data.individual_data.head()
```

### Find values with a methylation difference above a threshold

```python
DIFF_THRESHOLD = 0.01
large_diffs = data.diffs[(data.diffs >= DIFF_THRESHOLD).any(1)].dropna(how="all", axis=1)
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
        interesting[index] = row.loc[(row != 0)].dropna()
```

```python
import pandas as pd
```

```python
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
import plotly.express as px
```

### Raw quma results

```python
data.quma_results[data.means.index[0]].values.head()
```

```python
data.individual_data.head()
```

### Histograms of reads at a cell line and genomic position

```python
CELL = data.means.index[1]
POS = data.means.columns[5]
data.histogram(CELL, POS)
```

### Plotting read methylation distributions

```python
import pandas as pd
```

```python
df = data.individual_data.copy()
```

```python
df.head()
```

```python
def to_1D(series):
    """From https://towardsdatascience.com/dealing-with-list-values-in-pandas-dataframes-a177e534f173"""
    if isinstance(series, pd.Series):
        series = series.dropna()
        return pd.Series([x for _list in series for x in _list])
```

```python
POS = data.means.columns[5]
```

```python
to_1D(df.loc[:,POS]).value_counts()
```

```python
to_1D(df.loc[:,POS])
```

```python
to_1D(df.loc[:,POS]).mode()
```

#### Distribution at one cell line, position

```python
px.bar(to_1D(df.loc[:,POS]).value_counts(normalize=True))
```

#### Distribution across a cell line

```python
CELL = data.means.index[1]
```

```python
df2 = df.loc[CELL]
```

```python
df2.head()
```

```python
df2.explode().value_counts()
```

```python
px.violin(df2.explode())
```

#### Distribution across all cell lines

```python
all_values = []
```

```python
for each in df.index:
    temp = pd.Series(df.loc[each].explode())
#     print(temp)
    all_values.append(temp)
#     print(all_values)
```

```python
flat_list = [item for sublist in all_values for item in sublist]
```

```python
flat_series = pd.Series(flat_list)
```

```python
flat_series.value_counts()
```

```python
px.violin(flat_series)
```

### Identify large differences

```python
def find_big_diffs(df, min_diff: float):
    out = {}
    for each in df.columns:
    #     print(each)
        mean = to_1D(df[each].dropna()).mean()
        mode = to_1D(df[each].dropna()).mode()
    #     print(mode)
        if not mode.empty:   
    #         print(f"Mean: {mean}, mode: {mode}")
            diff = abs(mean - mode.values[0])
            if diff > min_diff:
                out[each] = diff
#                 print(f"Position: {each}, diff: {diff}")
#                 print(out)
    if out:
        return out
```

```python
big_ones = find_big_diffs(df, 0.1)
big_ones
```

## Statistical Tests for Normality

```python
summ = data.summarize_allelic_data()
```

```python
summ.head()
```

```python
summ.pivot(index="position", columns="cellLine", values="ad_stat")
```
