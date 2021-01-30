#!/usr/bin/env python3
"""pyllelic: a tool for detection of allelic-specific varation in DNA sequencing.

## Example usage in ipython / jupyter notebook:

    import pyllelic

    pyllelic.set_up_env_variables(
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
    )

    pyllelic.main()  # runs every step all at once

----------------------------------

## Example exploratory / step-by-step use in ipython / jupyter notebook:

    import pyllelic

    pyllelic.set_up_env_variables(
        base_path="/Users/abonham/documents/test_allelic/",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
    )

    files_set = pyllelic.make_list_of_bam_files()  # finds bam files

    positions = pyllelic.index_and_fetch(files_set)  # index and creates bam_output folders/files

    pyllelic.genome_parsing()  # writes out genome strings in bam_output folders

    cell_types = extract_cell_types(files_set)  # pulls out the cell types available for analysis

    df_list = run_quma_and_compile_list_of_df(cell_types, filename)  # run quma, get dfs

    means_df = process_means(df_list, positions, files_set)  # process means data from dataframes

    modes_df = process_modes(df_list, positions, cell_types)  # process modes data from dataframes

    diff_df = find_diffs(means_df, modes_df)  # find difference between mean and mode

    write_means_modes_diffs(means_df, modes_df, diffs_df, filename)  # write output to excel files
"""

# Imports
import pandas as pd
import numpy as np
import pysam
import os
import sys
from skbio.alignment import StripedSmithWaterman
import plotly.graph_objects as go
import subprocess
from pathlib import Path
from scipy import stats
from tqdm.notebook import tqdm
from typing import List, Dict, Set, Optional
from . import config


def set_up_env_variables(
    base_path: str, prom_file: str, prom_start: str, prom_end: str, chrom: str
) -> None:
    """Helper method to set up all our environmental variables, such as for testing.

    Args:
        base_path (str): directory where all processing will occur, put .bam files
                         in "test" sub-directory in this folder
        prom_file (str): filename of genmic sequence of promoter region of interest
        prom_start (str): start position to analyze in promoter region
        prom_end (str): final position to analyze in promoter region
        chrom (str): chromosome promoter is located on
    """

    config.base_directory = Path(base_path)
    config.promoter_file = Path(base_path) / prom_file
    config.results_directory = Path(base_path) / "results"
    config.bam_directory = Path(base_path) / "bam_output"
    config.analysis_directory = Path(base_path) / "test"
    config.promoter_start = prom_start
    config.promoter_end = prom_end
    config.chromosome = chrom


##################################################################################
##################################################################################
def main() -> None:
    """Run a given set of Pyllelic analysis, using values from
    supplied environmental variables.
    """

    if sys.argv[1]:
        filename: str = sys.argv[1]
    else:
        filename = "output.xlsx"

    files_set: List[str] = make_list_of_bam_files()
    positions: List[str] = index_and_fetch(files_set)
    genome_parsing()
    cell_types: List[str] = extract_cell_types(files_set)
    df_list: Dict[str, pd.DataFrame] = run_quma_and_compile_list_of_df(
        cell_types, filename
    )
    means_df: pd.DataFrame = process_means(df_list, positions, cell_types)
    modes_df: pd.DataFrame = process_modes(df_list, positions, cell_types)
    diffs_df: pd.DataFrame = find_diffs(means_df, modes_df)
    write_means_modes_diffs(means_df, modes_df, diffs_df, filename)


##################################################################################
##################################################################################


def genome_range(
    position: str, genome_string: str, offset: Optional[int] = None
) -> str:
    """Helper to return a genome string (e.g., "ATCGACTAG")
    given a position and an entire string.

    Args:
        position: genomic position on chromesome
        genome_string: string representation of genomic promoter known sequence
        offset: genomic position of promoter to offset reads

    Returns:
        str: genomic bases for indicated read / position
    """

    OFFSET: int = 1298163  # TERT offset

    if not offset:
        offset = OFFSET

    start: int = offset - (int(position) + 30)
    end: int = offset - (int(position) + 1)

    return genome_string[start:end]


def make_list_of_bam_files() -> List[str]:
    """Check analysis directory for all valid .bam files.

    Returns:
        list[str]: list of files
    """

    indv_bam: List[Path] = [b for b in config.analysis_directory.iterdir()]

    init_files_set: List[str] = [bam.name for bam in indv_bam if bam.suffix == ".bam"]

    bai_set: List[str] = [bai.name for bai in indv_bam if bai.suffix == ".bai"]

    f_set: List[str] = [str(bam).split(".")[0] for bam in init_files_set]

    baii: List[str] = [str(bai).split(".")[0] for bai in bai_set]

    files_set: List[str] = [name + (".TERT.bam") for name in f_set if name in baii]

    return files_set


def index_and_fetch(files_set: List[str]) -> List[str]:
    """Wrapper to call processing of each sam file.

    Args:
        files_set (list[str]): list of bam/sam files

    Returns:
        list[str]: list of genomic positions analyzed
    """

    sam_path: List[Path] = [config.base_directory / "test" / f for f in files_set]

    all_pos: Set = set()
    for sams in sam_path:
        pos: pd.Index = run_sam_and_extract_df(sams)
        all_pos.update(pos)

    return sorted(list(all_pos))


def run_sam_and_extract_df(sams: Path) -> pd.Index:
    """Process samfiles, pulling out sequence and position data
    and writing to folders/files.

    Args:
        sams (Path): path to a samfile

    Returns:
        pd.Index: list of unique positions in the samfile
    """

    # Index samfile if index file doesn't exist
    if not sams.with_suffix(".bai").exists():
        _: str = samtools_index(sams)  # we don't care what the output is

    # Grab the promoter region of interest
    samm: pysam.AlignmentFile = pysam.AlignmentFile(sams, "rb")
    itern = samm.fetch(
        config.chromosome, int(config.promoter_start), int(config.promoter_end)
    )

    position: List = []
    sequence: List = []

    for x in itern:
        cols = str(x).split()
        position.append(cols[3])
        sequence.append(cols[9])

    df: pd.DataFrame = pd.DataFrame(
        list(zip(position, sequence)), columns=["positions", "sequence"]
    )

    df2: pd.DataFrame = df.set_index("positions")
    # will set the inital index (on the leftmost column) to be position
    df3: pd.DataFrame = df2.stack()
    # if confused, see: https://www.w3resource.com/pandas/dataframe/dataframe-stack.php

    write_bam_output_files(sams, df2.index.unique(), df3)

    return df2.index.unique()


def write_bam_output_files(sams: Path, positions: List[str], df: pd.DataFrame) -> None:
    """Extract alignments from sequencing reads and output text files
       in bam directory.

    Args:
        sams (Path): path to a sam file
        positions (List[str]): list of unique positions
        df (pd.DataFrame): dataframe of sequencing reads
    """

    sam_name: str = sams.name

    # Make sure bam_output directory and sam subdirectories exist
    config.base_directory.joinpath("bam_output", sam_name).mkdir(
        parents=True, exist_ok=True
    )
    for each1 in positions:
        alignments: List = []

        # Set up query using alignment algorithm
        query_sequence: List[str] = df.loc[each1].head(1).tolist()[0]
        query: StripedSmithWaterman = StripedSmithWaterman(query_sequence)

        # Set up sequences to check for alignment
        target_sequences: List[str] = df.loc[each1].tolist()
        for target_sequence in target_sequences:
            alignment = query(target_sequence)
            alignments.append(alignment)

        read_file: List[str] = []
        for index, each in enumerate(alignments):
            read_file.append(str(">read" + str(index)))
            read_file.append(each.aligned_target_sequence)
            # returns aligned target sequence

        write_individual_bam_file(sam_name, each1, read_file)


def write_individual_bam_file(
    sam_name: str, filename: str, file_contents: List[str]
) -> None:
    """Write the contents of each bam output file.

    Args:
        sam_name (str): [description]
        filename (str): [description]
        file_contents (List[str]): [description]
    """
    directory: Path = config.base_directory.joinpath("bam_output", sam_name)
    with open(directory.joinpath(filename + ".txt"), "w") as file_handler:
        for item in file_contents:
            file_handler.write("{}\n".format(item))


def samtools_index(sams: Path) -> str:
    """Helper function to run external samtools index tool.

    Args:
        sams (Path): filepath to samfile

    Returns:
        str: output from samtools index shell command, usually discarded
    """

    command: List[str] = ["samtools", "index", os.fspath(sams)]

    output: subprocess.CompletedProcess = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: str = output.stdout

    return out


def genome_parsing(subfolders: List[Path] = None) -> None:
    """Writes out a list of genomic sequence strings for comparison to read data."""

    # Grab list of directories
    if not subfolders:
        subfolders = [x for x in config.bam_directory.iterdir() if x.is_dir()]

    # Grab genomic sequence
    with open(config.promoter_file, "r") as f:
        genome_base = f.readlines()
        genome_base_lines: List[str] = [s.rstrip("\n") for s in genome_base]
        genome_string: str = "".join(map(str, genome_base_lines))

    # Wrap everything in processing them one at a time
    for folder in tqdm(subfolders, desc="Cell Lines: "):

        # Grab list of read files in that directory:
        raw_read_files: List[str] = os.listdir(folder)
        read_files: List[str] = [
            os.path.splitext(i)[0]
            for i in raw_read_files
            if not i.lstrip().startswith("g")
            and not i.lstrip().startswith(".ip")
            and not i.lstrip().startswith(".DS")  # mac .DS_store files
        ]  # replace with regex?

        # Now, process each file:
        for read_name in read_files:
            file_lines: List = []
            # Grab the genomic sequence
            file_lines.append(str(">genome" + str(read_name)))
            file_lines.append(str(genome_range(read_name, genome_string)))

            # Save the reads as a text file for each position
            with open(
                folder.joinpath("g_" + str(read_name) + ".txt"), "w"
            ) as file_handler:
                for item in file_lines:
                    file_handler.write("{}\n".format(item))


def run_quma(directory: str, genomic_seq_file: str, reads_seq_file: str) -> str:
    """Helper function to run external QUMA tool.

    Args:
        directory (str): directory path to analyze
        genomic_seq_file (str): text file with known genomic sequence
        reads_seq_file (str): text file with experimental reads

    Returns:
        bytes: shell output from quma command
    """

    quma_path: str = os.fspath(config.base_directory.joinpath("quma_cui"))
    command: List[str] = [
        "perl",
        f"{quma_path}/quma.pl",
        "-g",
        f"{directory}/{genomic_seq_file}",
        "-q",
        f"{directory}/{reads_seq_file}",
    ]

    output: subprocess.CompletedProcess = subprocess.run(
        command, text=True, capture_output=True, check=True
    )
    out: str = output.stdout
    return out


def quma_full(cell_types, filename):
    """Run external QUMA methylation analysis on all specified cell lines.

    Args:
        cell_types (list[str]): list of cell lines in our dataset
        filename (str): desired output filename for xlsx output
    """

    # Grab list of directories
    subfolders: List[Path] = [p for p in config.bam_directory.iterdir() if p.is_dir()]

    writer: pd.ExcelWriter = pd.ExcelWriter(config.base_directory.joinpath(filename))

    # Wrap everything in processing them one at a time
    for folder in tqdm(subfolders, desc="Cell Lines"):

        # Get short name of cell_line
        cell_line_name: str = folder.name.split("_")[1]

        if set([cell_line_name]).intersection(set(cell_types)):
            # Set up a holding data frame from all the data
            holding_df: pd.DataFrame = pd.DataFrame()

            # Grab list of read files in that directory:
            raw_read_files: List[str] = os.listdir(folder)
            read_files: List[str] = [
                os.path.splitext(i)[0]
                for i in raw_read_files
                if i.endswith(".txt") and not i.lstrip().startswith("g")
            ]

            # Now, process each file:
            for read_name in tqdm(read_files, desc="Positions", leave=False):
                # file_lines = []

                quma_result: str = run_quma(
                    folder, "g_{}.txt".format(read_name), "{}.txt".format(read_name)
                )

                dots: List[str] = []
                for line in quma_result.splitlines():
                    if not line.lstrip().startswith("g"):
                        fields: List[str] = line.split("\t")
                        if float(fields[7]) < 80:
                            dots.append("FAIL")
                        else:
                            dots.append(fields[13])
                # Sort dots output by number of "1"
                dots2: List[str] = sorted(dots, key=lambda t: t.count("1"))

                # Next, add this readname to the holding data frame
                int_df: pd.DataFrame = pd.DataFrame({read_name: dots2})
                holding_df: pd.DataFrame = pd.concat([holding_df, int_df], axis=1)

            # Now, save it to an excel file
            holding_df.to_excel(writer, sheet_name=cell_line_name)

            del holding_df

    writer.save()


def extract_cell_types(file_sets: List[str]) -> List[str]:
    """Returns a list[str] of cell line names in the dataset."""

    return [file.split("_")[1] for file in file_sets]


def run_quma_and_compile_list_of_df(
    cell_types: List[str], filename: str, run_quma: bool = True
) -> Dict[str, pd.DataFrame]:
    """Wrapper to run QUMA on all cell lines in the dataset and write output file.

    Args:
        cell_types (list[str]): list of cell lines in the dataset
        filename (str): desired output filename
        run_quma (bool): whether we should invoke the external quma run, default True

    Returns:
        Dict[str, pd.DataFrame]: dataframe of quma results
    """

    if run_quma:
        quma_full(cell_types, filename)

    df: Dict[str, pd.DataFrame] = read_df_of_quma_results(filename)

    return df


def read_df_of_quma_results(filename: str) -> Dict[str, pd.DataFrame]:
    """Read excel output of quma results into dataframe.

    Args:
        filename (str): filename of quma results

    Returns:
        Dict[str, pd.DataFrame]: dict of dataframes of quma results
    """

    return pd.read_excel(
        config.base_directory.joinpath(filename),
        dtype=str,
        sheet_name=None,
        index_col=0,
    )


def process_means(
    dict_of_dfs: Dict[str, pd.DataFrame], positions: List[str], cell_types: List[str]
) -> pd.DataFrame:
    """Process the mean values at each position for each cell line.

    Args:
        dict_of_dfs (Dict[str, pd.DataFrame]): dict of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        pd.DataFrame: dataframes of mean values for each position in each cell line
    """

    # Gives the means of each positions-- NOT the mean of the entire dataframe!
    working_df: pd.DataFrame = pd.DataFrame()
    for pos in positions:
        working_df[pos] = ""
        for key, each in dict_of_dfs.items():
            values_list: List[float] = return_read_values(pos, key, dict_of_dfs)

            if values_list:
                pos_means: float = np.mean(values_list)
            else:  # No data or data doesn't meet minimums for analysis
                pos_means = np.nan

            working_df.loc[key, pos] = pos_means

    means_df: pd.DataFrame = working_df

    return means_df


def process_modes(
    dict_of_dfs: Dict[str, pd.DataFrame], positions: List[str], cell_types: List[str]
) -> pd.DataFrame:
    """Process the mode values at each position for each cell line.

    Args:
        dict_of_dfs (Dict[str, pd.DataFrame]): dict of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        pd.DataFrame: dataframes of mode values for each position in each cell line
    """

    # Gives the modes of each positions-- NOT the mean of the entire dataframe!
    working_df: pd.DataFrame = pd.DataFrame()
    for pos in positions:
        working_df[pos] = ""
        for key, each in dict_of_dfs.items():
            values_list: List[float] = return_read_values(pos, key, dict_of_dfs)

            if values_list:
                pos_modes: float = stats.mode(values_list)[0][0]
            else:  # No data or data doesn't meet minimums for analysis
                pos_modes = np.nan

            working_df.loc[key, pos] = pos_modes

    modes_df: pd.DataFrame = working_df

    return modes_df


def return_individual_data(
    dict_of_dfs: Dict[str, pd.DataFrame], positions: List[str], cell_types: List[str]
) -> pd.DataFrame:
    """Return a dataframe for methylation values at each position for each cell line.

    Args:
        dict_of_dfs (Dict[str, pd.DataFrame]): dict of dataframes of quma results
        positions (list[str]): list of genomic positions to analyze
        cell_types (list[str]): list of cell lines in the dataset

    Returns:
        pd.DataFrame: methylation values for each position in each cell line
    """

    working_df: pd.DataFrame = pd.DataFrame()
    for pos in tqdm(positions, desc="Position"):
        working_df[pos] = ""  # Create position column in dataframe
        for key, each in tqdm(dict_of_dfs.items(), desc="Cell Line", leave=False):
            values_list: List[float] = return_read_values(pos, key, dict_of_dfs)
            if values_list:
                data_for_df: List[float] = values_list
            else:  # No data or data doesn't meet minimums for analysis
                data_for_df = np.nan

            working_df.loc[key, pos] = data_for_df

    return working_df


def return_read_values(
    pos: str,
    key: str,
    dict_of_dfs: Dict[str, pd.DataFrame],
    min_reads: int = 4,
    min_sites: int = 2,
) -> List[float]:
    """Given a genomic position, cell line, and data, find fractional methylation
       values for each read in the dataset.

    Args:
        pos (str): genomic position of read
        key (str): cell line name
        dict_of_dfs (Dict[str, pd.DataFrame]): dictionary of methylation data
        min_reads (int): minimum bound for # reads to consider
        min_sites (int): minimum bound for # methylation sites to consider

    Returns:
        List[float]: list of fractional methylation for each read
    """

    bad_values: List[str] = ["N", "F"]  # for interpreting quma returns
    min_num_of_reads: int = min_reads  # Only include if >4 reads
    min_num_meth_sites: int = min_sites  # Only include is read has >2 methylation sites

    values_list: List[float] = []
    # Check if this cell line has that position
    if pos in dict_of_dfs.get(key).columns:  # type: ignore[union-attr]

        values_to_check: List[str] = (
            dict_of_dfs.get(key).loc[:, pos].dropna().astype(str)  # type: ignore[union-attr]
        )
        if not len(values_to_check) == 0:  # Skip empty values
            if (
                len(values_to_check) > min_num_of_reads
                and len(values_to_check) > min_num_meth_sites
            ):
                for value in values_to_check:
                    if not any(substring in value for substring in bad_values):
                        fraction_val: float = float(value.count("1")) / float(
                            len(value)
                        )  # number of methylated sites in each read
                        values_list.append(fraction_val)

    return values_list


def find_diffs(means_df: pd.DataFrame, modes_df: pd.DataFrame) -> pd.DataFrame:
    """Find the differences between means and modes for each cell line at each position
       in means and modes data.

    Args:
        means_df (pd.DataFrame): dataframe of means values
        modes_df (pd.DataFrame): dataframe of modes values

    Returns:
        pd.DataFrame: dataframe of difference values
    """

    return means_df.subtract(modes_df)


def write_means_modes_diffs(
    means_df: pd.DataFrame, modes_df: pd.DataFrame, diff_df: pd.DataFrame, filename: str
) -> None:
    """Wite out files of means, modes, and diffs for future analysis.

    Args:
        means_df (pd.DataFrame): dataframe of means values
        modes_df (pd.DataFrame): dataframe of modes values
        diff_df (pd.DataFrame): dataframe of diff values
        filename (str): desired root filename
    """

    means_df.to_excel(
        config.base_directory / str(filename + "_means.xlsx"),
        sheet_name="Means",
        index=True,
    )
    modes_df.to_excel(
        config.base_directory / str(filename + "_modes.xlsx"),
        sheet_name="Modes",
        index=True,
    )
    diff_df.to_excel(
        config.base_directory / str(filename + "_diff.xlsx"),
        sheet_name="Diff",
        index=True,
    )


def create_histogram(data: pd.DataFrame, cell_line: str, position: str) -> go.Figure:
    """Generate a graph figure showing fractional methylation in
       a given cell line at a given site.

    Args:
        data (pd.DataFrame): dataframe of individual data
        cell_line (str): name of cell line
        position (str): genomic position

    Returns:
        go.Figure: plotly figure object
    """
    fig: go.Figure = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=data.loc[cell_line, position],
            xbins=dict(
                start=-0.1,
                end=1.1,
                size=0.2,
            ),  # offset bins to center displayed bars
        )
    )
    fig.update_layout(
        title_text=f"Methylation at pos: {position} for cell line: {cell_line}",
        xaxis_title_text="Fraction of Sites Methylated in Read",
        yaxis_title_text="Read Count",
        bargap=0.2,
        template="seaborn",
    )

    return fig


def histogram(data: pd.DataFrame, cell_line: str, position: str) -> None:
    """Display a graph figure showing fractional methylation in
       a given cell line at a given site.

    Args:
        data (pd.DataFrame): dataframe of individual data
        cell_line (str): name of cell line
        position (str): genomic position
    """

    fig: go.Figure = create_histogram(data, cell_line, position)
    fig.show()


if __name__ == "__main__":
    main()
