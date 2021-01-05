#!/usr/bin/env python3
"""pyllelic: a tool for detection of allelic-specific varation in DNA sequencing.

## Example usage in ipython / jupyter notebook:

    import pyllelic

    pyllelic.set_up_env_variables(
        methyl_base="/Users/abonham/documents/methyl_test",
        promoter_seq="TERT-promoter-genomic-sequence.txt",
        promoter_start="1293000",
        promoter_end="1296000",
        chromosome="5",
    )

    pyllelic.main("output.xlsx")  # runs every step all at once

----------------------------------

## Example usage from command line:

    $ export METHYL_BASE=/Users/abonham/documents/methyl_test
    $ export PROMOTER_SEQ=TERT-promoter-genomic-sequence.txt
    $ export PROMOTER_START=1293000
    $ export PROMOTER_END=1296000
    $ export CHROMOSOME=5
    $ python pyllelic.py output.xlsx

----------------------------------

## Example exploratory / step-by-step use in ipython / jupyter notebook:

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

    cell_types = extract_cell_types(files_set)  # pulls out the cell types available for analysis

    df_list = run_quma_and_compile_list_of_df(cell_types, filename)  # run quma, get dfs

    means_df = process_means(df_list, positions, files_set)  # process means data from dataframes

    write_means_to_excel(means_df, files_set)  # write means data to excel files

"""

# Imports
import pandas as pd
import numpy as np
import pysam
import os
import sys
from skbio.alignment import StripedSmithWaterman
import plotly.express as px  # noqa
import subprocess
from pathlib import Path
from io import StringIO
import statistics as stat  # noqa
import openpyxl as pxl

# Define base constants, but run set_up_env_variables() to correctly update them at run time
BASE_DIRECTORY = (
    Path(os.environ.get("METHYL_BASE")) if os.environ.get("METHYL_BASE") else Path.cwd()
)
BASE_DIRECTORY.mkdir(exist_ok=True)

PROMOTER_FILE = (
    BASE_DIRECTORY / os.environ.get("PROMOTER_SEQ")
    if os.environ.get("PROMOTER_SEQ")
    else BASE_DIRECTORY / "promoter.txt"
)

RESULTS_DIRECTORY = BASE_DIRECTORY / "results"
RESULTS_DIRECTORY.mkdir(exist_ok=True)

BAM_DIRECTORY = BASE_DIRECTORY / "bam_output"
BAM_DIRECTORY.mkdir(exist_ok=True)

ANALYSIS_DIRECTORY = (
    Path(os.environ.get("ANALYSIS_FOLDER"))
    if os.environ.get("ANALYSIS_FOLDER")
    else BASE_DIRECTORY / "test"
)
ANALYSIS_DIRECTORY.mkdir(exist_ok=True)

PROMOTER_START = (
    int(os.environ.get("PROMOTER_START"))
    if os.environ.get("PROMOTER_START")
    else 1293000
)

PROMOTER_END = (
    int(os.environ.get("PROMOTER_END")) if os.environ.get("PROMOTER_END") else 1296000
)

CHROMOSOME = os.environ.get("CHROMOSOME") if os.environ.get("CHROMOSOME") else "5"


def set_up_env_variables(
    methyl_base="/Users/abonham/documents/methyl_test",
    promoter_seq="TERT-promoter-genomic-sequence.txt",
    promoter_start="1293000",
    promoter_end="1296000",
    chromosome="5",
):
    """Helper method to set up all our environmental variables, such as for testing."""

    os.environ["METHYL_BASE"] = methyl_base
    os.environ["PROMOTER_SEQ"] = promoter_seq
    os.environ["PROMOTER_START"] = promoter_start
    os.environ["PROMOTER_END"] = promoter_end
    os.environ["CHROMOSOME"] = chromosome

    BASE_DIRECTORY = (
        Path(os.environ.get("METHYL_BASE"))
        if os.environ.get("METHYL_BASE")
        else Path.cwd()
    )
    BASE_DIRECTORY.mkdir(exist_ok=True)

    PROMOTER_FILE = (  # noqa
        BASE_DIRECTORY / os.environ.get("PROMOTER_SEQ")
        if os.environ.get("PROMOTER_SEQ")
        else BASE_DIRECTORY / "promoter.txt"
    )

    RESULTS_DIRECTORY = BASE_DIRECTORY / "results"
    RESULTS_DIRECTORY.mkdir(exist_ok=True)

    BAM_DIRECTORY = BASE_DIRECTORY / "bam_output"
    BAM_DIRECTORY.mkdir(exist_ok=True)

    ANALYSIS_DIRECTORY = (
        Path(os.environ.get("ANALYSIS_FOLDER"))
        if os.environ.get("ANALYSIS_FOLDER")
        else BASE_DIRECTORY / "test"
    )
    ANALYSIS_DIRECTORY.mkdir(exist_ok=True)

    PROMOTER_START = (  # noqa
        int(os.environ.get("PROMOTER_START"))
        if os.environ.get("PROMOTER_START")
        else 1293000
    )

    PROMOTER_END = (  # noqa
        int(os.environ.get("PROMOTER_END"))
        if os.environ.get("PROMOTER_END")
        else 1296000
    )

    CHROMOSOME = os.environ.get("CHROMOSOME") if os.environ.get("CHROMOSOME") else "5"  # noqa


##################################################################################
##################################################################################
def main(filename):
    """Run a given set of Pyllelic analysis, using values from supplied environmental variables.

    Args:
        filename: filename to write output of analysis to.
    """

    files_set = make_list_of_bam_files()
    positions = index_and_fetch(files_set)
    genome_parsing()
    cell_types = extract_cell_types(files_set)
    df_list = run_quma_and_compile_list_of_df(cell_types, filename)
    means_df = process_means(df_list, positions, cell_types)
    write_means_to_excel(means_df, cell_types)


##################################################################################
##################################################################################


def genome_range(position, genome_string):
    """Helper to return a genome string (e.g., "ATCGACTAG")
    given a position and an entire string.

    Args:
        position: genomic position on chromesome
        genome_string: string representation of genomic promoter known sequence

    Returns:
        str: genomic bases for indicated read / position
    """

    start = 1298163 - (int(position) + 30)
    end = 1298163 - (int(position) + 1)

    return genome_string[start:end]


def run_quma(directory, genomic_seq_file, reads_seq_file):
    """Helper function to run external QUMA tool.

    Args:
        directory (str): directory path to analyze
        genomic_seq_file (str): text file with known genomic sequence
        reads_seq_file (str): text file with experimental reads

    Returns:
        bytes: shell output from quma command
    """

    quma_path = os.fspath(BASE_DIRECTORY.joinpath("quma_cui"))
    command = [
        "perl",
        f"{quma_path}/quma.pl",
        "-g",
        f"{directory}/{genomic_seq_file}",
        "-q",
        f"{directory}/{reads_seq_file}",
    ]

    out = subprocess.run(command, text=True, capture_output=True).stdout
    return out


def make_list_of_bam_files():
    """Check analysis directory for all valid .bam files.

    Returns:
        list[str]: list of files
    """

    indv_bam = [b for b in ANALYSIS_DIRECTORY.iterdir()]

    init_files_set = [bam.name for bam in indv_bam if bam.suffix == ".bam"]

    bai_set = [bai.name for bai in indv_bam if bai.suffix == ".bai"]

    f_set = [str(bam).split(".")[0] for bam in init_files_set]

    baii = [str(bai).split(".")[0] for bai in bai_set]

    files_set = [file + (".TERT.bam") for file in f_set if file in baii]

    final_file_sets = []

    # Code truncates to aleviate excel file_name limit of <=31 characters:
    final_file_sets = [data[:29] for data in files_set]

    return final_file_sets


def index_and_fetch(files_set):
    """Wrapper to call processing of each sam file.

    Args:
        files_set (list[str]): list of bam/sam files

    Returns:
        list[str]: list of genomic positions analyzed
    """

    sam_path = [BASE_DIRECTORY / "test" / f for f in files_set]

    all_pos = []
    for sams in sam_path:
        pos = run_sam_and_extract_df(sams)
        all_pos.append(pos)

    return sorted(all_pos)


def run_sam_and_extract_df(sams):
    """Process samfiles, pulling out sequence and position data and writing to folders/files.

    Args:
        sams (str): path to a samfile

    Returns:
        list: list of unique positions in the samfile
    """

    # Make sure each sam file has an index by calling external samtools index function
    _ = samtools_index(sams)  # we don't care what the output is

    # Grab the promoter region of interest
    samm = pysam.AlignmentFile(sams, "rb")
    itern = samm.fetch(CHROMOSOME, PROMOTER_START, PROMOTER_END)

    position = []
    sequence = []

    for x in itern:
        cols = str(x).split()
        position.append(cols[3])
        sequence.append(cols[9])

    # Transfer into dataframe for processing
    df = pd.DataFrame(list(zip(position, sequence)), columns=["positions", "sequence"])

    df2 = df.set_index("positions")
    # will set the inital index (on the leftmost column) to be position
    df3 = df2.stack()
    # if confused, see: https://www.w3resource.com/pandas/dataframe/dataframe-stack.php

    # Now, iterate through the dataframe
    for each1 in df2.index.unique():
        alignments = []

        # Set up query using alignment algorithm
        query_sequence = df3.loc[each1].head(1).tolist()[0]
        query = StripedSmithWaterman(query_sequence)

        # Set up sequences to check for alignment
        target_sequences = df3.loc[each1].tolist()
        for target_sequence in target_sequences:
            alignment = query(target_sequence)
            alignments.append(alignment)

        read_file = []
        for index, each in enumerate(alignments):
            read_file.append(str(">read" + str(index)))
            read_file.append(alignments[0].aligned_target_sequence)
            # returns aligned target sequence

        # Make sure bam_output directory and sam subdirectories exist
        BASE_DIRECTORY.joinpath("bam_output", sams.name).mkdir(
            parents=True, exist_ok=True
        )
        directory = BASE_DIRECTORY.joinpath("bam_output", sams.name)

        with open(directory.joinpath(str(each1) + ".txt"), "w") as file_handler:
            for item in read_file:
                file_handler.write("{}\n".format(item))

    return df2.index.unique()  # return all unique positions in the data


def samtools_index(sams):
    """Helper function to run external samtools index tool.

    Args:
        sams (str): filepath to samfile

    Returns:
        str: output from samtools index shell command, usually discarded
    """

    command = ["samtools", "index", os.fspath(sams)]

    out = subprocess.run(command, capture_output=True, text=True)  # remove shell=True,
    return out


def genome_parsing():
    """Writes out a list of genomic sequence strings for comparison to read data."""

    # Grab list of directories
    subfolders = [x for x in BAM_DIRECTORY.iterdir() if x.is_dir()]

    # Grab genomic sequence
    with open(PROMOTER_FILE, "r") as f:
        genome_base = f.readlines()
        genome_base_lines = [s.rstrip("\n") for s in genome_base]
        genome_string = "".join(map(str, genome_base_lines))

    # Wrap everything in processing them one at a time
    for folder in subfolders:

        # Grab list of read files in that directory:
        raw_read_files = os.listdir(folder)
        read_files = [
            os.path.splitext(i)[0]
            for i in raw_read_files
            if not i.lstrip().startswith("g")
            and not i.lstrip().startswith(".ip")
            and not i.lstrip().startswith(".DS")
        ]

        # Now, process each file:
        for read_name in read_files:
            file_lines = []
            # Grab the genomic sequence and write it
            file_lines.append(str(">genome" + str(read_name)))
            file_lines.append(str(genome_range(read_name, genome_string)))

            # Save the reads as a text file for each position
            with open(
                folder.joinpath("g_" + str(read_name) + ".txt"), "w"
            ) as file_handler:
                for item in file_lines:
                    file_handler.write("{}\n".format(item))


def quma_full(cell_types, filename):
    """Run external QUMA methylation analysis on all specified cell lines.

    Args:
        cell_types (list[str]): list of cell lines in our dataset
        filename (str): desired output filename for xlsx output
    """

    # Grab list of directories
    subfolders = [f.path for f in os.scandir(BAM_DIRECTORY) if f.is_dir()]

    writer = pd.ExcelWriter(BASE_DIRECTORY.joinpath(filename))

    # Wrap everything in processing them one at a time
    for folder in subfolders:

        if any(substring in folder for substring in cell_types):
            # Set up a holding data frame from all the data
            holding_df = pd.DataFrame()

            # Grab list of read files in that directory:
            raw_read_files = os.listdir(folder)
            read_files = [
                os.path.splitext(i)[0]
                for i in raw_read_files
                if i.endswith(".txt") and not i.lstrip().startswith("g")
            ]

            # Now, process each file:
            for read_name in read_files:
                # file_lines = []

                quma_result = run_quma(
                    folder, "g_{}.txt".format(read_name), "{}.txt".format(read_name)
                )

                dots = []
                for line in StringIO(quma_result.decode("utf-8")):
                    if not line.lstrip().startswith("g"):
                        fields = line.split("\t")
                        if float(fields[7]) < 80:
                            dots.append("FAIL")
                        else:
                            dots.append(fields[13])
                # Sort dots output by number of "1"
                dots2 = sorted(dots, key=lambda t: t.count("1"))

                # Next, add this readname to the holding data frame
                int_df = pd.DataFrame({read_name: dots2})
                holding_df = pd.concat([holding_df, int_df], axis=1)
                # holdng_df[read_name]=dots

            # Now, save it to an excel file
            holding_df.to_excel(writer, sheet_name=Path(folder).name)

            del holding_df

    writer.save()


def extract_cell_types(file_sets):
    """Returns a list[str] of cell lines in the dataset."""

    return [file.split("_")[1] for file in file_sets]


def run_quma_and_compile_list_of_df(cell_types, filename):
    """Wrapper to run QUMA on all cell lines in the dataset and write output files.

    Args:
        cell_types (list[str]): list of cell lines in the dataset
        filename (str): desired output filename

    Returns:
        list[pd.DataFrame]: list of dataframes of quma results
    """

    df_full_list = []

    quma_full(
        cell_types, filename
    )
    df = pd.read_excel(BASE_DIRECTORY.joinpath(filename), dtype=str, sheet_name=None)
    df_full_list.append(df)

    for df in df_full_list:
        for each in df:
            df[each].replace(to_replace="nan", value=np.nan, inplace=True)

    return df_full_list


def process_means(list_of_dfs, positions, cell_types):
    """Process the mean values at each position for each cell line.

    Args:
        list_of_dfs (list[pd.DataFrame]): list of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        list[pd.DataFrame]): list of dataframes of mean values for each position in each cell line
    """

    bad_values = ["N", "F"]  # for interpreting quma returns

    df_full = list_of_dfs

    total_means = []

    # Gives the means of each individual positions-- NOT the mean of the entire dataframe!

    means_cols = []
    # FIXME: Need to specify keys correctly, should be dictionary instead?
    means_cols.append(list(df_full[0].keys()))
    for pos in positions:
        mean_col = []
        mean_col.append(pos)
        for ix, each in enumerate(df_full):
            values_list = []
            if pos in df_full[each].columns:
                if not (
                    len(df_full[each].loc[:, pos].dropna().astype(str)) < 5
                    and not len(df_full[each].loc[:, pos].dropna().astype(str)[0]) < 3
                ):

                    for value in df_full[each].loc[:, pos].dropna().astype(str):
                        if not any(substring in value for substring in bad_values):
                            values_list.append(
                                float(value.count("1")) / float(len(value))
                            )
                        else:
                            pass
                else:
                    pass
            else:
                pass
            a = np.mean(values_list)  # TODO: THIS IS THE DIFFERENT LINE
            try:
                b = a.item()
            except ValueError:
                b = np.nan
            mean_col.append(b)

        means_cols.append(mean_col)

        means_df = pd.DataFrame(means_cols).transpose()

        alpha = "cell_types"  # FIXME: not actually unique naming

        means_df.to_csv(path_or_buf=RESULTS_DIRECTORY.joinpath(alpha + pos))
        total_means.append(means_df)

    return total_means


def write_means_to_excel(means_df, cell_types):
    """Write means data to a compiled xlsx file.

    Args:
        means_df (list[pd.DataFramme]): list of means dataframes for all the cell lines
        cell_types (list[str]): list of genomic positions analyzed
    """

    for ix, each in enumerate(means_df):
        t = cell_types[ix]  # FIXME: how to access correct cell_type label??
        sheet = t
        excel_book = pxl.load_workbook(RESULTS_DIRECTORY.joinpath("testT.xlsx"))
        # FIXME: Fixed filename??
        with pd.ExcelWriter("testT.xlsx", engine="openpyxl") as writer:
            writer.book = excel_book
            writer.sheets = {
                worksheet.title: worksheet for worksheet in excel_book.worksheets
            }

            means_df.to_excel(writer, sheet, index=True)
        writer.save()


def new_process_means():
    pass


"""Remaining Code:

nums=-1
df_full_total = []
for df in df_full:
    num = nums + 1
    file_name = file_sets[num]
    z = df_full[file_name]
    df_full_total.append(z)
   # df_full['fh_CALU1_LUNG.TERT.bam'][1:]


---------

indiv_mean = []
n = -1
num_list = []
Mean_name = ['',]
Mean_data = []
for h in means_df.transpose():
    n = n + 1
    num_list.append(n)


for n in num_list:
    means_d = means_df.loc[n]
    indiv_mean.append(means_d)
for sets in indiv_mean:
    sets = sets[0]
    sets=str(sets)
    if sets == str(sets):
        #print(sets)
        Mean_name.append(sets)
#dic_Mean = pd.DataFrame(indiv_mean)
# for sets in indiv_mean:
for sets in indiv_mean:

    sets = sets[1:]
    Mean_data.append(sets)
#print(Mean_name)
#print(Means_name)
#print(Mean_name)
#print(Mean_data)
#print(MeaN)
MeaN = dict(zip(Mean_name,Mean_data))

Means_name = pd.DataFrame(MeaN.items())

# Mean_data = pd.DataFrame(list(Mean_data.items()))
# Means_position

means_excel = pd.DataFrame(MeaN).transpose()
means_excel
# dic_Mean = df.set_index('id')
# dic_Mean.transpose()

# combined = pd.merge(Means_name, Mean_data, left_index=True, right_index=True)

------------

total_mean = []
total_mode = []
total_percent_diff = []



means_np = means_df.to_numpy()[1:]

means_lis = means_np.tolist()
k = []
y = -1
for mean in means_lis:
    mean=mean[1:]

    mean1 =[]

    for h in mean:
        if h == float(h):

            mean1.append(h)

    final_mode = stat.mode(mean1)
    final_mean = stat.mean(mean1)
    total_mean.append(final_mean)
    total_mode.append(final_mode)


    avg = ((final_mode + final_mean)/2)
#avg = ((final_mode + final_mean)/2)
    per_diff = abs((final_mode - final_mean)*100/avg)
    total_percent_diff.append(per_diff)

# print(total_mode)
# print(total_mean)
#print(total_percent_diff)

-----------

diff= []

for a in total_percent_diff:
    diff.append(a)

dic_mean = dict(zip(files_set,total_mean),orient ='index')
dic_mode = dict(zip(files_set,total_mode),orient ='index')

dif_ex = dict(zip(files_set,diff),orient ='index')

diff = pd.DataFrame(list(dif_ex.items()))
dic_mean = pd.DataFrame(list(dic_mean.items()))
dic_mode = pd.DataFrame(list(dic_mode.items()))


combined = pd.merge(dic_mode, dic_mean, left_index=True, right_index=True)
all_comb = pd.merge(diff, combined, left_index=True, right_index=True)

all_comb

--------

a = 'diff_' +files_set[0]
sheet = a
excel_book = pxl.load_workbook(base_path.joinpath('testT.xlsx'))
with pd.ExcelWriter('testT.xlsx', engine='openpyxl') as writer:
    writer.book = excel_book
    writer.sheets = {
        worksheet.title: worksheet
        for worksheet in excel_book.worksheets
    }

    secondMocksDF = all_comb
    secondMocksDF.to_excel(writer, sheet, index=True)

writer.save()

------------

## Modes Histogram
counts, bins = np.histogram(modes_df.T[1:].set_index(0).values, bins=np.linspace(0,1,5))
bins = 0.5 * (bins[:-1] + bins[1:])
fig_modes = px.bar(x=bins, y=counts, labels={'x':'Methylation Percent', 'y':'count'})
fig_modes.show()

newdf = means_df
keepdf =newdf.rename(columns=newdf.iloc[0]).drop(newdf.index[0]) #renames columns.
keepdf = keepdf.set_index('')
keepdf = keepdf.transpose().dropna()
keepdf = pd.DataFrame(keepdf)
print(keepdf)
binned_data = []
for col in keepdf.columns:
    col_bin_means = stats.binned_statistic(
                                       col,
                                       keepdf[col].values.tolist(),
                                       bins=5,
                                       range=(0, 1)
     )[0]
     a = np.histogram(keepdf[col].values.tolist(), bins=5, range=(0,1), density=True)[0]
#creates histogram of iterative objects in keepdf
     a_True = list(a / a.sum())
#gives a list of percent of each object by dividing each value by the sum of the values.
    a_True.insert(0,str(col))
    binned_data.append(a_True)
binned_data

mut_binned_df = (
    pd.DataFrame(binned_data,columns=['position','0%','25%','50%','75%','100%'])
    .set_index('position')
    )
mut_binned = mut_binned_df.dropna().transpose()
mut_binned
mut_binned_df.values.tolist()
"""

if __name__ == "__main__":
    main(sys.argv[1])  # run the whole shebang using defaults and env variables!
