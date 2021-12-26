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

<img src="https://github.com/Paradoxdruid/pyllelic/blob/master/assets/pyllelic_logo.png?raw=true" alt="pyllelic logo" width="100"/>


## Background


This notebook illustrates the import and use of `pyllelic` in a jupyter environment.

Source code: <https://github.com/Paradoxdruid/pyllelic>

Documentation: <https://paradoxdruid.github.io/pyllelic/>


## Pre-setup / File preparation


### Obtaining fastq data


We can download rrbs (reduced representation bisulfite sequencing) data from the Encode project:
<http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/>


Those files are in unaligned fastq format.  We will need to align these to a reference human genome.


### Retrieve promoter sequence


We will download the genomic region of interest from the UCSC Genome browser and save it in a text file:

<!-- #region -->
```python
from pyllelic import process
# Retrieve promoter genomic sequence
process.retrieve_promoter_seq("{prom_filename}.txt", chrom: "chr5", start: 1293200, end: 1296000)
```
<!-- #endregion -->

### Preparing bisultife-converted genome


To prepare the genome and align reads, we'll use [bismark](https://github.com/FelixKrueger/Bismark), bowtie2, and samtools (through its pysam wrapper).

First, we need to download a genomic index sequence: <http://hgdownload.soe.ucsc.edu/goldenPath/hg19>

<!-- #region -->
```python
# Processing imports
from pathlib import Path

# Set up file paths
genome = Path("/{your_directory}/{genome_file_directory}")
fastq = Path("/{your_directory}/{your_fastq_file.fastq.gz}")
```
<!-- #endregion -->

**WARNING:** The next command is processor, RAM, and time intensive, and only needs to be run once!

<!-- #region -->
```python
# Prepare genome via bismark
process.prepare_genome(genome) # can optionally give path to bowtie2 if not in PATH
```
<!-- #endregion -->

### Aligning reads


**WARNING:** The next command is processor, RAM, and time intensive, and only needs to be run once!

<!-- #region -->
```python
# Convert fastq to bismark-aligned bam
from pyllelic import process
process.bismark(genome, fastq)
```
<!-- #endregion -->

Notes:

* We recommend renaming your `.bam` file to encode cell-line and tissue.  Our convention is: `fh_CELLLINE_TISSUE.REGION.bam`

Next, we need to sort and index the bam file using samtools functions.

<!-- #region -->
```python
# Sort the bamfile
bamfile = Path("/{your_directory}/{bam_filename}.bam")
process.pysam_sort(bamfile)
```
<!-- #endregion -->

<!-- #region -->
```python
# Create an index of the sorted bamfile
sorted_bam = Path("/{your_directory}/{bam_filename}_sorted.bam")
process.pysam_index(sorted_bam)
```
<!-- #endregion -->

### Organize directory contents


* Place sorted bam files and index files (again, rename to capture cell-line and tissue info) in the `test` folder for analysis by pyllelic.
* Place the promoter sequence in your main directory.


## Set-up

```python
from pyllelic import pyllelic
```

1. Set up your disk location:  ```base_path``` should be the directory we'll do our work in
2. Make a sub-directory under ```base_path``` with a folder named ```test``` and put the ```.bam``` and ```.bai``` files in ```test```

```python
config = pyllelic.configure(
    base_path="/home/andrew/allellic/",  # Windows WSL set-up
    #     base_path="/Users/abonham/documents/test_allelic/",  # OSX setup
    prom_file="tert_genome.txt",
    prom_start="1293200",
    prom_end="1296000",
    chrom="5",
    offset=1293000,
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
files_set = files_set[0:10]
```

### Perform full methylation analysis and generate data object

**Warning**: This step is processor and time intensive and will perform all processing and analysis of the bam file data.

```python
data = pyllelic.GenomicPositionData(config=config, files_set=files_set)
```

### Check main data outcomes

```python
", ".join(data.positions)
```

```python
data.means.head()
```

```python
data.modes.head()
```

```python
data.diffs.head()
```

```python
data.histogram("JHOM1", "1295089")
```

```python
data.diffs.loc[:,data.diffs.any()]
```

```python
data.quma_results[data.means.index[0]].values.head()
```

```python
data.individual_data.loc[data.means.index[0], data.means.columns[0]]
```

## Write Output

### Save entire object as pickle

```python
# data.save_pickle("test_data.pickle")
```

### Save Quma Results to excel

```python
# data.save("output.xlsx")
```

### Save analysis files (means, modes, diffs) to excel

```python
# data.write_means_modes_diffs("Full_") # Sets the filename stub
```

### Reopen saved object

```python
# data = pyllelic.GenomicPositionData.from_pickle("test_data2.pickle")
```

## Data Analysis

### View raw data of methylation percentage per read

```python
df = data.individual_data
```

### Statistical Tests for Normality

```python
summ = data.summarize_allelic_data()
```

```python
summ.head()
```

```python
ad_big_diffs = summ.pivot(index="position", columns="cellLine", values="ad_stat")
```

```python
ad_big_diffs = ad_big_diffs.dropna(axis=1, how="all").count(axis=1).to_frame()
```

```python
ad_big_diffs.head()
```

### Raw quma results

```python
data.quma_results[data.means.index[0]].values.head()
```

```python
data.individual_data.head()
```

## Visualizing Data


### Histograms of reads at a cell line and genomic position

```python
# CELL = data.means.index[13]
# POS = data.means.columns[15]
# data.histogram(CELL, POS)
```

### Heatmap of mean methylation values

```python
data.heatmap(min_values=40, width=800, height=2000, data_type="diffs")
```

### Bar chart of significant methylation differences

```python
data.sig_methylation_differences()
```

### Find values with a methylation difference above a threshold

```python
# import plotly.express as px
```

```python
import pandas as pd
```

```python
# DIFF_THRESHOLD = 0.1
```

```python
# large_diffs = data.diffs[(data.diffs >= DIFF_THRESHOLD).any(1)].dropna(
#     how="all", axis=1
# )
```

```python
# large_diffs.head()
```

```python
# interesting = {}
# for index, row in large_diffs.iterrows():
#     if row.loc[(row != 0)].dropna().any():
#         interesting[index] = row.loc[(row != 0)].dropna()
```

```python
# big_diffs = pd.DataFrame.from_dict(interesting)
```

```python
# big_diffs.head()
```

```python
# big_diffs_counts = big_diffs.dropna(axis=1, how="all").count(axis=1).to_frame()
```

```python
# fig = px.bar(big_diffs_counts, template="seaborn")
# fig.update_layout(xaxis_type="linear", showlegend=False, barmode="stack")
# fig.update_xaxes(
#     tickformat="r", tickangle=45, nticks=40, title="Position", range=[1292981, 1295979]
# )
# fig.update_yaxes(title="# of large differences")
# fig.update_traces(width=50)
# fig.show()
```

### Plotting read methylation distributions

```python
df = data.individual_data.copy()
```

```python
# df.head()
```

```python
def to_1D(series):
    """From https://towardsdatascience.com/dealing-with-list-values-in-pandas-dataframes-a177e534f173"""
    if isinstance(series, pd.Series):
        series = series.dropna()
        return pd.Series([x for _list in series for x in _list])
```

```python
POS = "1295089"  # data.means.columns[5]
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
# big_ones
```

```python
big_ones
```

```python
ad_big_diffs.shape
```

```python
# import plotly.express as px

# fig = px.bar(big_ones, template="seaborn")
# fig.update_layout(
#     xaxis_type="linear",
#     showlegend=False,
#     title="TERT Dataset Significant Methylation Differences",
#     barmode="overlay",
# )
# fig.update_xaxes(
#     tickformat="r", tickangle=45, nticks=40, title="Position", range=[1292981, 1295979]
# )
# fig.update_yaxes(title="# of significant differences")
# fig.update_traces(width=50)
# fig.show()
```

```python
ad_big_diffs.values.flatten()
```

```python
import plotly.graph_objects as go
```

```python
fig = go.Figure()
fig.add_trace(go.Bar(x=ad_big_diffs.index, y=ad_big_diffs.values.flatten()))
fig.update_layout(
    xaxis_type="linear",
    showlegend=False,
    title="Significant Methylation Differences",
#     barmode="overlay",
    template="seaborn"
)
fig.update_xaxes(
    tickformat="r", tickangle=45, nticks=40, title="Position", range=[int(ad_big_diffs.index.min()), int(ad_big_diffs.index.max())],
)
fig.update_yaxes(title="# of significant differences")
fig.update_traces(width=50)
fig.show()
```

## New graphs

```python
df = data.means.copy()
```

```python
import itertools
```

```python
N = len(df.index)
N
```

```python
df['grp'] = list(itertools.chain.from_iterable([x]*5 for x in range(0, int(N/5))))
```

```python
df.mean().head()
```

```python
df.groupby('grp').mean()
```

## New methylation graphs

```python
%load_ext nb_black
```

```python
df = data.individual_data
```

```python
df.head(1)
```

```python
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
```

```python
def make_binary(data):
    if isinstance(data, list):
        new = [data.count(0), data.count(1)]
    else:
        new = [0, 0]
    return new
```

```python
def make_methyl_df(df, row):
    new_df = pd.DataFrame(
        {each[0]: each[1].values[0] for each in df.loc[row].to_frame().iterrows()}
    )
    new_df = new_df.rename(index={0: "unmethylated", 1: "methylated"})
    return new_df
```

```python
def make_methyl_bar(df):
    fig = px.bar(
        df.T,
        height=20,
        width=800,
        color_discrete_sequence=["white", "black"],
        labels={"value": "reads", "variable": "status"},
    )
    fig.update_layout(showlegend=False, margin=dict(l=0, r=0, t=0, b=0))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig
```

```python
def make_it_all(df):
    df2 = df.applymap(make_binary)
    fig_dict = {}
    for each in df2.index:
        methyl_df = make_methyl_df(df2, each)
        new_fig = make_methyl_bar(methyl_df)
        fig_dict[each] = new_fig
    return fig_dict
```

```python
newest = make_it_all(df)
```

```python
import plotly.subplots as sp
```

```python
def make_stacked_fig(fig_dict):
    this_figure = sp.make_subplots(
        rows=len(fig_dict),
        cols=2,
        column_widths=[0.1, 0.9],
        horizontal_spacing=0,
        vertical_spacing=0.005,
    )

    for i, (key, each) in enumerate(fig_dict.items()):
        this_figure.append_trace(each["data"][0], row=i + 1, col=2)
        this_figure.append_trace(each["data"][1], row=i + 1, col=2)
        this_figure.add_annotation(
            text=key,
            xref="x",
            yref="y",
            x=1,
            y=1,
            showarrow=False,
            font=dict(
                size=12,
            ),
            row=i + 1,
            col=1,
        )

    this_figure.update_layout(
        showlegend=False,
        margin=dict(l=0, r=0, t=40, b=0),
        barmode="stack",
        width=800,
        height=len(fig_dict) * 40,
        template="ggplot2",
        paper_bgcolor="rgba(132,132,132,1)",
        title_text="Methylated and Unmethylated Reads by Cell Line",
        title={"font": {"color": "white"}},
    )
    this_figure.update_xaxes(visible=False)
    this_figure.update_yaxes(visible=False)

    return this_figure
```

```python
full_stack = make_stacked_fig(newest)
```

```python
full_stack.show()
```

```python

```

```python

```
