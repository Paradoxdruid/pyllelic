# pyllelic

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Paradoxdruid/pyllelic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Paradoxdruid/pyllelic/context:python)  [![CodeFactor](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic/badge)](https://www.codefactor.io/repository/github/paradoxdruid/pyllelic)  [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

**pyllelic**: a tool for detection of allelic-specific varation in DNA sequencing.

:exclamation: **This is a work-in-progress, and may not be functional at present!** :exclamation:

## Example usage in ipython / jupyter notebook:
```python
    import pyllelic

    pyllelic.set_up_env_variables(
        methyl_base="/Users/abonham/documents/methyl_test",
        promoter_seq="TERT-promoter-genomic-sequence.txt",
        promoter_start="1293000",
        promoter_end="1296000",
        chromosome="5",
    )

    pyllelic.main("output.xlsx")  # runs every step all at once
```

----------------------------------

## Example usage from command line:

```bash
    $ export METHYL_BASE=/Users/abonham/documents/methyl_test
    $ export PROMOTER_SEQ=TERT-promoter-genomic-sequence.txt
    $ export PROMOTER_START=1293000
    $ export PROMOTER_END=1296000
    $ export CHROMOSOME=5
    $ python pyllelic.py output.xlsx
```

----------------------------------

## Example exploratory / step-by-step use in ipython / jupyter notebook:

```python
    import pyllelic

    pyllelic.set_up_env_variables(
        methyl_base="/Users/abonham/documents/methyl_test",
        promoter_seq="TERT-promoter-genomic-sequence.txt",
        promoter_start="1293000",
        promoter_end="1296000",
        chromosome="5",
    )

    files_set = pyllelic.make_list_of_bam_files()  # finds bam files

    positions = pyllelic.index_and_fetch(files_set)  # index bam and creates bam_output folders/files

    pyllelic.genome_parsing()  # writes out genome strings in bam_output folders

    cell_types = pyllelic.extract_cell_types(files_set)  # pulls out the cell types available for analysis

    df_list = pyllelic.run_quma_and_compile_list_of_df(cell_types, filename)  # run quma, get dfs

    means_df = pyllelic.process_means(df_list, positions, files_set)  # process means data from dataframes

    pyllelic.write_means_to_excel(means_df, files_set)  # write means data to excel files
```
