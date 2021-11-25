# pyllelic

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Paradoxdruid/pyllelic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Paradoxdruid/pyllelic/context:python)  [![CodeFactor](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic/badge)](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic)  [![Codacy Badge](https://app.codacy.com/project/badge/Grade/c8c86fe25a644cb69b8b6e789ca1c18f)](https://www.codacy.com/gh/Paradoxdruid/pyllelic/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Paradoxdruid/pyllelic&amp;utm_campaign=Badge_Grade)  [![Codacy Badge](https://app.codacy.com/project/badge/Coverage/c8c86fe25a644cb69b8b6e789ca1c18f)](https://www.codacy.com/gh/Paradoxdruid/pyllelic/dashboard)  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)  [![PyPI](https://img.shields.io/pypi/v/pyllelic?color=success)](https://pypi.org/project/pyllelic/) ![GitHub](https://img.shields.io/github/license/Paradoxdruid/fealden)

<p align="right">
  ⭐ &nbsp;&nbsp;the project to show your appreciation. :arrow_upper_right:
</p>

<img src="https://github.com/Paradoxdruid/pyllelic/blob/master/assets/pyllelic_logo.png" width="100" height="100" style="float: left; margin-right: 10px;">

**pyllelic**: a tool for detection of allelic-specific methylation variation in bisulfite DNA sequencing files.

Pyllelic documention is available at **<https://paradoxdruid.github.io/pyllelic/>** and see [`pyllelic_notebook.ipynb`](https://github.com/Paradoxdruid/pyllelic/blob/master/pyllelic_notebook.md) for an interactive demonstration.

## Example exploratory use in jupyter notebook

```python
    from pyllelic import pyllelic

    config = pyllelic.configure(  # Specify file and directory locations
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293200",
        prom_end="1296000",
        chrom="5",
    )

    files_set = pyllelic.make_list_of_bam_files(config)  # finds bam files

    # Run pyllelic; make take some time depending on number of bam files
    data = pyllelic.GenomicPositionData(config=config, files_set=files_set)

    positions = data.positions

    cell_types = data.cell_types

    means_df = data.means

    modes_df = data.modes
    
    diff_df = data.diffs

    individual_data = data.individual_data

    data.save("output.xlsx")  # save methylation results

    data.save_pickle("my_run.pickle")  # save data object for later analysis
    
    data.write_means_modes_diffs(filename="Run1_")  # write output data files

    data.histogram("CELL_LINE", "POSITION")  # visualize data for a point

    data.heatmap(min_values=1)  # methylation level heatmap

    data.quma_results["CELL_LINE"]  # see summary data for a cell line
```

----------------------------------

## Dependencies and Installation

### Conda environment

Create a new conda environment using python 3.7:

```bash
conda create --name PYLLELIC python=3.7
conda activate PYLLELIC
```

### Install pyllelic

```bash
pip install pyllelic
```

or

```bash
git clone https://github.com/Paradoxdruid/pyllelic.git
```

## Authors

This software is developed as academic software by [Dr. Andrew J. Bonham](https://github.com/Paradoxdruid) at the [Metropolitan State University of Denver](https://www.msudenver.edu). It is licensed under the GPL v3.0.

This software incorporates implementation from [QUMA](http://quma.cdb.riken.jp), licensed under the GPL v3.0.
