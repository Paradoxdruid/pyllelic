# pyllelic

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Paradoxdruid/pyllelic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Paradoxdruid/pyllelic/context:python)  [![CodeFactor](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic/badge)](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic)  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

:microscope: **pyllelic**: a tool for detection of allelic-specific methylation variation in DNA sequencing files.

:warning: This is a work-in-progress, and not all functions are implemented at present. :warning:

See [`pyllelic_notebook.ipynb`](https://github.com/Paradoxdruid/pyllelic/blob/master/pyllelic_notebook.ipynb) for an interactive demonstration.

## Example usage in ipython / jupyter notebook:
```python
    import pyllelic

    pyllelic.set_up_env_variables(
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
    )

    pyllelic.main("output.xlsx")  # runs every step all at once
```

----------------------------------

## Example exploratory / step-by-step use in ipython / jupyter notebook:

```python
    import pyllelic

    pyllelic.set_up_env_variables(  # Specify file and directory locations
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
    )

    pyllelic.setup_directories()  # Read env variables to set up directories to use

    files_set = pyllelic.make_list_of_bam_files()  # finds bam files

    positions = pyllelic.index_and_fetch(files_set)  # index bam and creates bam_output folders/files

    pyllelic.genome_parsing()  # writes out genome strings in bam_output folders

    cell_types = pyllelic.extract_cell_types(files_set)  # pulls out the cell types available for analysis

    df_list = pyllelic.run_quma_and_compile_list_of_df(cell_types, filename)  # run quma, get dfs

    means_df = pyllelic.process_means(df_list, positions, files_set)  # process means data from dataframes

    modes_df = pyllelic.process_modes(df_list, positions, cell_types)  # process modes data from dataframes
    
    diff_df = pyllelic.find_diffs(means_df, modes_df)  # find difference between mean and mode

    pyllelic.write_means_modes_diffs(means_df, modes_df, diffs_df, filename)  # write output data to excel files
```

----------------------------------

## Dependencies and Installation
### Conda Environment
* Create a new conda environment using python 3.7:
```bash
conda create --name methyl python=3.7
conda activate methyl
conda config --add channels conda-forge
conda config --set channel_priority strict
```
* Install python dependencies:
```bash
conda install pandas numpy scipy plotly dash notebook xlsxwriter xlrd
conda install -c bioconda samtools pysam scikit-bio
```
* Install system dependencies:
```bash
conda install -c bioconda emboss
conda install -c bioconda perl perl-app-cpanminus
cpan install Statistics::Lite
```
* Set up jupyter:
```bash
conda install -c conda-forge jupyter_contrib_nbextensions
```
### Install quma
* Download from http://quma.cdb.riken.jp/files/quma_cui-1.0.0.tar.gz
* Untar into `base_path`


## Authors
This software is developed as academic software by [Dr. Andrew J. Bonham](https://github.com/Paradoxdruid) at the [Metropolitan State University of Denver](https://www.msudenver.edu). It is licensed under the GPL v3.0.
