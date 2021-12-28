# pyllelic

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Paradoxdruid/pyllelic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Paradoxdruid/pyllelic/context:python)  [![CodeFactor](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic/badge)](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic)  [![Codacy Badge](https://app.codacy.com/project/badge/Grade/c8c86fe25a644cb69b8b6e789ca1c18f)](https://www.codacy.com/gh/Paradoxdruid/pyllelic/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Paradoxdruid/pyllelic&amp;utm_campaign=Badge_Grade)  [![Codacy Badge](https://app.codacy.com/project/badge/Coverage/c8c86fe25a644cb69b8b6e789ca1c18f)](https://www.codacy.com/gh/Paradoxdruid/pyllelic/dashboard)  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)  

[![PyPI](https://img.shields.io/pypi/v/pyllelic?color=success)](https://pypi.org/project/pyllelic/) [![Anaconda-Server Badge](https://anaconda.org/paradoxdruid/pyllelic/badges/version.svg)](https://anaconda.org/paradoxdruid/pyllelic) ![GitHub](https://img.shields.io/github/license/Paradoxdruid/fealden)

<img align="left" src="https://github.com/Paradoxdruid/pyllelic/blob/master/assets/pyllelic_logo.png?raw=true" width="100" height="100">

**pyllelic**: a tool for detection of allelic-specific methylation variation in bisulfite DNA sequencing files.

Pyllelic documention is available at **<https://paradoxdruid.github.io/pyllelic/>** and see [`pyllelic_notebook.ipynb`](https://github.com/Paradoxdruid/pyllelic/blob/master/pyllelic_notebook.md) for an interactive demonstration.

## Quickstart

Run an interactive sample pyllelic environment in your web browser using [`mybinder.org`](https://mybinder.org):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Paradoxdruid/pyllelic/HEAD?filepath=pyllelic_binder_example.ipynb)

## pyllelic in action

<p align="center">
<img src="https://github.com/Paradoxdruid/pyllelic/blob/master/assets/pyllelic_demo.gif?raw=true" alt="pyllelic demo gif">
</p>

## Dependencies and Installation

### Using Conda (preferred)

Create a new conda environment using python 3.8:

Easiest:

```bash
# Get environment.yml file from this repo
curl -L https://github.com/Paradoxdruid/pyllelic/blob/master/environment.yml?raw=true > env.yml

# Create and activate conda environment
conda env create --file=env.yml
conda activate pyllelic
```

<details>
  <summary>or more explictly step by step instructions</summary>

```bash
conda create --name pyllelic python=3.8
conda activate pyllelic
conda config --env --add channels conda-forge
conda config --env --add channels bioconda
conda config --env --add channels paradoxdruid
conda install pyllelic 

# Optional but usual use case:
conda install notebook jupyter_contrib_nbextensions ipywidgets
```

</details>

### Docker container

```bash
docker pull ghcr.io/paradoxdruid/pyllelic:latest
```

### PyPi installation

<details>
  <summary>PyPi instructions</summary>

This will require independent installation of samtools, bowtie2, and bismark packages.

```bash
# PyPi
python3 -m pip install pyllelic
# or Github
python3 -m pip install git+https://github.com/Paradoxdruid/pyllelic.git
```

</details>

## Example exploratory use in jupyter notebook

```python
    from pyllelic import pyllelic

    config = pyllelic.configure(  # Specify file and directory locations
        base_path="/home/jovyan/assets/",
        prom_file="tert_genome.txt",
        prom_start="1293200",
        prom_end="1296000",
        chrom="5",
        offset=1293000, # start position of retrieved promoter sequence
    )

    files_set = pyllelic.make_list_of_bam_files(config)  # finds bam files

    # Run pyllelic; make take some time depending on number of bam files
    data = pyllelic.GenomicPositionData(config=config, files_set=files_set)

    positions = data.positions

    cell_types = data.cell_types

    means_df = data.means  # mean methylation of reads

    modes_df = data.modes  # mode methylation of reads
    
    diff_df = data.diffs  # difference mean - mode of reads

    individual_data = data.individual_data  # read methylation values

    data.save("output.xlsx")  # save methylation results

    data.save_pickle("my_run.pickle")  # save data object for later analysis
    
    data.write_means_modes_diffs(filename="Run1_")  # write output data files

    data.histogram("CELL_LINE", "POSITION")  # visualize data for a point

    data.heatmap(min_values=1)  # methylation level heatmap

    data.quma_results["CELL_LINE"]  # see summary data for a cell line
```

----------------------------------

## Authors

This software is developed as academic software by [Dr. Andrew J. Bonham](https://github.com/Paradoxdruid) at the [Metropolitan State University of Denver](https://www.msudenver.edu). It is licensed under the GPL v3.0.

This software incorporates implementation from [QUMA](http://quma.cdb.riken.jp), licensed under the GPL v3.0.
