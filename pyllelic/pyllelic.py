#!/usr/bin/env python3
"""pyllelic: a tool for detection of allelic-specific variation
   in reduced representation bisulfate DNA sequencing.
"""

import os
import signal
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, NamedTuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pysam
from Bio import pairwise2
from scipy import stats
from tqdm.notebook import tqdm

from . import quma
from .config import Config

import logging

logging.basicConfig(filename="pyllelic_test.log", level=logging.DEBUG)

# Initialize shared configuration object
config = Config()

# Initialized multiprocessing limits
NUM_THREADS = cpu_count() - 1


class AD_stats(NamedTuple):
    """Helper class for NamedTuple results from anderson_darling_test"""

    sig: bool
    stat: float
    crits: List[Any]


class BamOutput:
    """Storage container to process BAM sequencing files and store processed results."""

    def __init__(self, sam_directory: Path, genome_string: str, config: Config) -> None:
        self._config: Config = config
        self.name: str = str(sam_directory)
        self.values: Dict[str, str] = {}
        self.positions: pd.Index = self.run_sam_and_extract_df(sam_directory)
        self.genome_values: Dict[str, str] = {}
        self.genome_parsing(genome_string)

    def run_sam_and_extract_df(self, sams: Path) -> pd.Index:
        """Process samfiles, pulling out sequence and position data
        and writing to folders/files.

        Args:
            sams (Path): path to a samfile

        Returns:
            pd.Index: list of unique positions in the samfile
        """

        # Index samfile if index file doesn't exist
        if not sams.with_suffix(".bai").exists():
            _: bool = self.pysam_index(sams)  # we don't care what the output is

        # Grab the promoter region of interest
        samm: pysam.AlignmentFile = pysam.AlignmentFile(sams, "rb")
        itern = samm.fetch(
            self._config.chromosome,
            int(self._config.promoter_start),
            int(self._config.promoter_end),
        )

        position: List[str] = []
        sequence: List[str] = []

        for x in itern:
            cols = str(x).split()
            position.append(cols[3])
            sequence.append(cols[9])

        df: pd.DataFrame = pd.DataFrame(
            list(zip(position, sequence)), columns=["positions", "sequence"]
        )

        df2: pd.DataFrame = df.set_index("positions")
        df3: pd.DataFrame = df2.stack()

        self.write_bam_output(df2.index.unique(), df3)

        return df2.index.unique()

    def write_bam_output(self, positions: List[str], df: pd.DataFrame) -> None:
        """Extract alignments from sequencing reads and output text strings
        in bam values dictionary.

        Args:
            sams (Path): path to a sam file
            positions (List[str]): list of unique positions
            df (pd.DataFrame): dataframe of sequencing reads
        """

        for each1 in positions:
            alignments: List[str] = []

            # Set up query using alignment algorithm
            query_sequence: List[str] = df.loc[each1].head(1).tolist()[0]

            # Set up sequences to check for alignment
            target_sequences: List[str] = df.loc[each1].tolist()
            for target_sequence in target_sequences:
                alignment = pairwise2.align.localms(
                    query_sequence, target_sequence, 2, -3, -5, -2
                )
                aligned_segment = alignment[0].seqA[
                    alignment[0].start : alignment[0].end
                ]
                alignments.append(aligned_segment)

            read_file: List[str] = []
            for index, each in enumerate(alignments):
                read_file.append(str(">read" + str(index)))
                read_file.append(each)
                # returns aligned target sequence

            self.values[each1] = "\n".join(read_file)

    @staticmethod
    def pysam_index(bamfile: Path) -> bool:
        """Helper function to run external samtools index.

        Args:
            bamfile (Path): filepath to bam file

        Returns:
            bool: verification of samtools command, usually discarded
        """

        pysam.index(os.fspath(bamfile))

        return True

    @staticmethod
    def _genome_range(
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

        # OFFSET: int = 1298163  # TERT offset

        if not offset:
            offset = config.offset

        start: int = offset - (int(position) + 30)
        end: int = offset - (int(position) + 1)

        return genome_string[start:end]

    def genome_parsing(self, genome_string: str) -> None:
        """Writes out a list of genomic sequence strings for comparison to read data."""

        # Grab list of read files in that directory:
        read_files: List[str] = list(self.values.keys())

        # Now, process each file:
        for read_name in read_files:
            file_lines: List[str] = []
            # Grab the genomic sequence
            file_lines.append(str(">genome" + str(read_name)))
            file_lines.append(str(self._genome_range(read_name, genome_string)))

            # Save the reads as a text file for each position
            self.genome_values[read_name] = "\n".join(file_lines)


class QumaResult:
    """Storage container to process and store quma-style methylation results."""

    def __init__(
        self, read_files: List[str], genomic_files: List[str], positions: List[str]
    ) -> None:
        self._read_files: List[str] = read_files
        self._genomic_files: List[str] = genomic_files
        self._positions: List[str] = positions
        self.values: pd.DataFrame = self._pool_processing()

    @staticmethod
    def process_raw_quma(quma_result: str) -> List[str]:
        """Spit and process raw quma results into a pandas dataframe.

        Args:
            quma_result (str): console output from quma program.

        Returns:
            pd.DataFrame: intermediate dataframe to append
        """

        dots: List[str] = []
        for line in quma_result.splitlines():
            if not line.lstrip().startswith("g"):
                fields: List[str] = line.split("\t")
                if float(fields[7]) < 70:  # checks for high mercent match
                    dots.append("FAIL")
                else:
                    dots.append(fields[13])
        # Sort dots output by number of "1"
        return sorted(dots, key=lambda t: t.count("1"))

    def _pool_processing(self) -> pd.DataFrame:
        """Helper to run quma in a multiprocessing pool.

        Returns:
            pd.DataFrame: dataframe of quma results
        """

        # Set up a holding data frame from all the data
        holding_df: pd.DataFrame = pd.DataFrame()

        returns: List[Any] = []
        with Pool(NUM_THREADS, self._init_worker) as pool:

            for position, read, genomic in zip(
                self._positions, self._read_files, self._genomic_files
            ):

                result = pool.apply_async(
                    self._thread_worker,
                    (genomic, read, position),
                )
                returns.append(result)

            results = [result.get() for result in returns]

        # Add to excel file
        for each in results:
            try:
                holding_df = pd.concat([holding_df, each], axis=1)
            except (AttributeError, KeyError):
                pass

        return holding_df

    @staticmethod
    def _init_worker() -> None:
        """
        Pool worker initializer for keyboard interrupt on Windows
        """
        signal.signal(signal.SIGINT, signal.SIG_IGN)

    def _thread_worker(
        self, genomic_contents: str, read_contents: str, position: str
    ) -> pd.DataFrame:
        """Queue worker for quma functions.

        Args:

        Returns:
            pd.DataFrame: dataframe of quma results
        """

        quma_result: str = self.access_quma(genomic_contents, read_contents)

        processed_quma: List[str] = self.process_raw_quma(quma_result)
        # Next, add this readname to the holding data frame
        int_df: pd.DataFrame = pd.DataFrame({position: processed_quma})
        return int_df

    @staticmethod
    def access_quma(genomic_contents: str, read_contents: str) -> str:
        """Helper function to run internal QUMA tool.

        Args:

        Returns:
            str: output from quma command
        """

        result = quma.quma_main(genomic_contents, read_contents)

        return result


##################################################################################


class GenomicPositionData:
    """Class to process reduced representation bisulfite methylation sequencing data.

    When initialized, GenomicPositionData reads sequencing file (.bam) locations from
    a config object, and then automatically performs alignment into BamOutput objects,
    and then performs methylation analysis, storing the results as QumaResults.

    Finally, the aggregate data is analyzed to create some aggregate metrics such as
    means, modes, and differences (diffs), as well as expose methods for plotting
    and statistical analysis.
    """

    def __init__(self, config: Config, files_set: List[str]) -> None:
        self.config: Config = config
        self.files_set: List[str] = files_set

        with open(self.config.promoter_file, "r") as f:
            genome_base = f.readlines()
            genome_base_lines: List[str] = [s.rstrip("\n") for s in genome_base]
            self.genome_string: str = "".join(map(str, genome_base_lines))

        self._bam_output: Dict[str, BamOutput] = {}
        self.index_and_fetch()
        self.positions: List[str] = self._calculate_positions()
        self.cell_types = self._bam_output.keys()

        self.quma_results: Dict[str, QumaResult] = self.quma_full_threaded()

        self.means: pd.DataFrame = self.process_means()
        self.modes: pd.DataFrame = self.process_modes()
        self.diffs: pd.DataFrame = self.find_diffs(self.means, self.modes)
        self.individual_data: pd.DataFrame = self.return_individual_data()

    def index_and_fetch(self) -> None:
        """Wrapper to call processing of each sam file.

        Args:
            genome_string (str): genomic sequence to align against

        Returns:
            list[str]: list of genomic positions analyzed
        """

        sam_path: List[Path] = [
            self.config.base_directory / "test" / f for f in self.files_set
        ]

        for sams in tqdm(sam_path):
            self._bam_output[sams] = BamOutput(sams, self.genome_string, self.config)

    def _calculate_positions(self) -> List[str]:
        """Return sorted list of all positions analyzed.

        Returns:
            List[str]: sorted list of all positions analyzed.
        """

        positions = []
        for sublist in self._bam_output.values():
            for each in sublist.values.keys():
                positions.append(each)

        return sorted(set(positions))

    @staticmethod
    def extract_cell_types(file_sets: List[str]) -> List[str]:
        """Returns a list[str] of cell line names in the dataset."""

        return [file.split("_")[1] for file in file_sets]

    def quma_full_threaded(self) -> Dict[str, QumaResult]:
        """Run external QUMA methylation analysis on all specified cell lines,
        using multiprocessing library.

        Returns:
            Dict[str, QumaResult]: dictionary of cell_types and quma results
        """

        quma_results: Dict[str, QumaResult] = {}
        for name, bam_result in tqdm(self._bam_output.items(), desc="Cell Lines"):
            cell_line_name: str = name.name.split("_")[1]

            read_files: List[str] = [each for each in bam_result.values.values()]
            genomic_files: List[str] = [
                each for each in bam_result.genome_values.values()
            ]
            positions = self.positions

            quma_results[cell_line_name] = QumaResult(
                read_files, genomic_files, positions
            )

        return quma_results

    def save(self, filename: str = "output.xlsx") -> None:
        """Save quma results to an excel file.

        Args:
            filename (str, optional): Filename to save to. Defaults to "output.xlsx".
        """

        writer: pd.ExcelWriter = pd.ExcelWriter(
            config.base_directory.joinpath(filename)
        )

        for name, each in self.quma_results.items():
            each.values.to_excel(writer, sheet_name=name)

        writer.save()

    def process_means(self) -> pd.DataFrame:
        """Process the mean values at each position for each cell line.

        Returns:
            pd.DataFrame: dataframes of mean values for each position in each cell line
        """

        # Gives the means of each positions-- NOT the mean of the entire dataframe!
        working_df: pd.DataFrame = pd.DataFrame()
        for pos in self.positions:
            working_df[pos] = ""
            for key in self.quma_results.keys():
                values_list: List[float] = self.return_read_values(pos, key)
                if values_list:
                    pos_means: float = float(np.mean(values_list))
                else:  # No data or data doesn't meet minimums for analysis
                    pos_means = np.nan

                working_df.at[key, pos] = pos_means

        means_df: pd.DataFrame = working_df

        return means_df

    def process_modes(self) -> pd.DataFrame:
        """Process the mode values at each position for each cell line.

        Returns:
            pd.DataFrame: dataframes of mode values for each position in each cell line
        """

        # Gives the modes of each positions-- NOT the mean of the entire dataframe!
        working_df: pd.DataFrame = pd.DataFrame()
        for pos in self.positions:
            working_df[pos] = ""
            for key in self.quma_results.keys():
                values_list: List[float] = self.return_read_values(pos, key)

                if values_list:
                    pos_modes: float = stats.mode(values_list)[0][0]
                else:  # No data or data doesn't meet minimums for analysis
                    pos_modes = np.nan

                working_df.at[key, pos] = pos_modes

        modes_df: pd.DataFrame = working_df

        return modes_df

    def return_individual_data(self) -> pd.DataFrame:
        """Return a dataframe for methylation values at each position for each cell line.

        Returns:
            pd.DataFrame: methylation values for each position in each cell line
        """

        working_df: pd.DataFrame = pd.DataFrame()
        for pos in tqdm(self.positions, desc="Position"):
            working_df[pos] = ""  # Create position column in dataframe
            for key in tqdm(self.quma_results.keys(), desc="Cell Line", leave=False):
                values_list: List[float] = self.return_read_values(pos, key)
                if values_list:
                    data_for_df: Union[List[float], float] = values_list
                else:  # No data or data doesn't meet minimums for analysis
                    data_for_df = np.nan

                working_df.at[key, pos] = data_for_df

        return working_df

    def return_read_values(
        self,
        pos: str,
        key: str,
        min_sites: int = 1,
    ) -> List[float]:
        """Given a genomic position, cell line, and data, find fractional methylation
        values for each read in the dataset.

        Args:
            pos (str): genomic position of read
            key (str): cell line name
            min_sites (int): minimum bound for # methylation sites to consider

        Returns:
            List[float]: list of fractional methylation for each read
        """

        bad_values: List[str] = ["N", "F"]  # for interpreting quma returns
        min_num_meth_sites: int = (
            min_sites  # Only include is read has >2 methylation sites
        )

        values_list: List[float] = []
        # Check if this cell line has that position
        df: pd.DataFrame = self.quma_results[key].values
        if pos in df.columns:  # type: ignore[union-attr]

            values_to_check: List[str] = self.get_str_values(df, pos)
            for value in values_to_check:
                if not any(substring in value for substring in bad_values):
                    if len(value) > min_num_meth_sites:
                        fraction_val: float = float(value.count("1")) / float(
                            len(value)
                        )  # number of methylated sites in each read
                        values_list.append(fraction_val)

        return values_list

    @staticmethod
    def get_str_values(df: pd.DataFrame, pos: str) -> List[str]:
        """Return list of values a a position in given dataframe of methylation data.

        Args:
            df (pd.DataFrame): cell line dataframe
            pos (str): position being analyzed

        Returns:
            List[str]: list of values
        """
        return df.loc[:, pos].dropna().astype(str)  # type: ignore[union-attr]

    @staticmethod
    def find_diffs(means_df: pd.DataFrame, modes_df: pd.DataFrame) -> pd.DataFrame:
        """Find the differences between means and modes for each cell line at each pos
        in means and modes data.

        Args:
            means_df (pd.DataFrame): dataframe of means values
            modes_df (pd.DataFrame): dataframe of modes values

        Returns:
            pd.DataFrame: dataframe of difference values
        """

        return means_df.subtract(modes_df)

    @staticmethod
    def truncate_diffs(diffs_df: pd.DataFrame) -> pd.DataFrame:
        """Remove missing or non-interesting diffs, and sort by magnitude.

        Args:
            diffs_df (pd.DataFrame): dataframe of diffs values

        Returns:
            pd.DataFrame: truncated dataframe
        """
        return diffs_df.dropna(how="all")

    def write_means_modes_diffs(
        self,
        filename: str,
    ) -> None:
        """Wite out files of means, modes, and diffs for future analysis.

        Args:
            filename (str): desired root filename
        """

        self.means.to_excel(
            config.base_directory / str(filename + "_means.xlsx"),
            sheet_name="Means",
            index=True,
        )
        self.modes.to_excel(
            config.base_directory / str(filename + "_modes.xlsx"),
            sheet_name="Modes",
            index=True,
        )
        self.diffs.to_excel(
            config.base_directory / str(filename + "_diff.xlsx"),
            sheet_name="Diff",
            index=True,
        )

    @staticmethod
    def create_histogram(
        data: pd.DataFrame, cell_line: str, position: str
    ) -> go.Figure:
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

        fig.update_xaxes(range=[-0.1, 1.1])

        return fig

    def histogram(self, cell_line: str, position: str) -> None:
        """Display a graph figure showing fractional methylation in
        a given cell line at a given site.

        Args:
            cell_line (str): name of cell line
            position (str): genomic position
        """
        data = self.individual_data

        fig: go.Figure = self.create_histogram(data, cell_line, position)
        fig.show()

    def generate_ad_stats(self) -> pd.DataFrame:
        """Generate Anderson-Darling normality statistics for an individual data df.

        Returns:
            pd.DataFrame: df of a-d test statistics
        """
        np.seterr(divide="ignore", invalid="ignore")  # ignore divide-by-zero errors
        df = self.individual_data
        df2 = df.applymap(self.anderson_darling_test)
        return df2

    def summarize_allelic_data(self) -> pd.DataFrame:
        """Create a dataframe only of likely allelic methylation positions

        Returns:
            pd.DataFrame: dataframe of cell lines with likely allelic positions
        """
        np.seterr(divide="ignore", invalid="ignore")  # ignore divide-by-zero errors
        sig_dict: Dict[str, List[Any]] = {
            "cellLine": [],
            "position": [],
            "ad_stat": [],
            "p_crit": [],
            "diff": [],
            "raw": [],
        }
        for index, row in self.individual_data.iterrows():
            for column in row.index:
                value = row[column]
                if np.all(pd.notnull(value)):
                    good, stat, crits = self.anderson_darling_test(value)
                    if good:
                        sig_dict["diff"].append(self.diffs.loc[index, column])
                        sig_dict["cellLine"].append(index)
                        sig_dict["position"].append(column)
                        sig_dict["raw"].append(value)
                        sig_dict["ad_stat"].append(stat)
                        sig_dict["p_crit"].append(crits[4])
        sig_df = pd.DataFrame(sig_dict)
        return sig_df

    @staticmethod
    def anderson_darling_test(raw_list: Optional[pd.Series]) -> AD_stats:
        """Run the Anderson-Darling normality test on methylation data at a point.

        Args:
            raw_list (pd.Series): list of fractional methylation levels per read.

        Returns:
            AD_stats[bool, float, List[Any]]: is the data significantly allelic (bool),
                                    A-D statistic (float), critical values (list)
        """

        if np.all(pd.notnull(raw_list)):
            stat: float
            crits: List[Any]
            stat, crits, _ = stats.anderson(raw_list)
            is_sig: bool = bool(stat > crits[4])
            return AD_stats(is_sig, stat, crits)
        return AD_stats(False, np.nan, [np.nan])


##################################################################################
##################################################################################

# Module level helper functions


def set_up_env_variables(
    base_path: str,
    prom_file: str,
    prom_start: str,
    prom_end: str,
    chrom: str,
    offset: int,
) -> None:
    """Helper method to set up all our environmental variables, such as for testing.

    Args:
        base_path (str): directory where all processing will occur, put .bam files
                         in "test" sub-directory in this folder
        prom_file (str): filename of genmic sequence of promoter region of interest
        prom_start (str): start position to analyze in promoter region
        prom_end (str): final position to analyze in promoter region
        chrom (str): chromosome promoter is located on
        offset (int): genomic position of promoter to offset reads
    """

    config.base_directory = Path(base_path)
    config.promoter_file = Path(base_path) / prom_file
    config.results_directory = Path(base_path) / "results"
    config.bam_directory = Path(base_path) / "bam_output"
    config.analysis_directory = Path(base_path) / "test"
    config.promoter_start = prom_start
    config.promoter_end = prom_end
    config.chromosome = chrom
    config.offset = offset


def make_list_of_bam_files() -> List[str]:
    """Check analysis directory for all valid .bam files.

    Returns:
        list[str]: list of files
    """
    return [f.name for f in config.analysis_directory.iterdir() if f.suffix == ".bam"]


# Main magic
if __name__ == "__main__":
    if sys.argv[1]:
        filename: str = sys.argv[1]
    else:  # pragma: no cover
        filename = "output.xlsx"

    files_to_analyze = make_list_of_bam_files()

    data = GenomicPositionData(config=config, files_set=files_to_analyze)
