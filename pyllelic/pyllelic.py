#!/usr/bin/env python3
"""pyllelic: a tool for detection of allelic-specific variation
   in reduced representation bisulfate DNA sequencing.
"""

import pickle  # nosec
import signal
from multiprocessing import Pool, cpu_count
from multiprocessing.pool import AsyncResult
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, TypeVar, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pysam
from Bio import pairwise2
from scipy import stats
from tqdm.auto import tqdm

from . import quma
from .config import Config

GPD = TypeVar("GPD", bound="GenomicPositionData")

# Initialized multiprocessing limits
NUM_THREADS = cpu_count() - 1

# SANITY ALIGNMENT LIMITS
MIN_NUMBER_READS = 1
MIN_PERCENT_ALIGNMENT = 50


class AD_stats(NamedTuple):
    """Helper class for NamedTuple results from anderson_darling_test"""

    sig: bool
    stat: float
    crits: List[np.ndarray]


class BamOutput:
    """Storage container to process BAM sequencing files and store processed results."""

    def __init__(self, sam_directory: Path, genome_string: str, config: Config) -> None:
        self._config: Config = config
        self.name: str = str(sam_directory)
        """str: path to bam file analyzed."""

        self.values: Dict[str, str] = {}
        """Dict[str, str]: dictionary of reads at a given position"""

        self.positions: pd.Index = self._run_sam_and_extract_df(sam_directory)
        """"pd.Index: index of genomic positions in the bam file."""

        self.genome_values: Dict[str, str] = {}
        """Dict[str, str]: dictionary of read files and contents."""

        self._genome_parsing(genome_string)

    def _run_sam_and_extract_df(self, sams: Path) -> pd.Index:
        """Process samfiles, pulling out sequence and position data
        and writing to folders/files.

        Args:
            sams (Path): path to a samfile

        Returns:
            pd.Index: list of unique positions in the samfile
        """

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
        df3: pd.Series = df2.stack()

        self._write_bam_output(df2.index.unique(), df3)

        return df2.index.unique().tolist()

    def _write_bam_output(self, positions: List[str], df: pd.Series) -> None:
        """Extract alignments from sequencing reads and output text strings
        in bam values dictionary.

        Args:
            sams (Path): path to a sam file
            positions (List[str]): list of unique positions
            df (pd.Series): series of sequencing reads
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

    def _genome_range(
        self, position: str, genome_string: str, offset: Optional[int] = None
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
            offset = self._config.offset

        start: int = offset - (int(position) + 30)
        end: int = offset - (int(position) + 1)
        if genome_string[start:end]:
            return genome_string[start:end]
        raise ValueError("Position not in genomic promoter file")

    def _genome_parsing(self, genome_string: str) -> None:
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
        """pd.DataFrame: dataframe of quma methylation analysis values."""

    @staticmethod
    def _process_raw_quma(quma_result: str) -> List[str]:
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
                if float(fields[7]) < MIN_PERCENT_ALIGNMENT:
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

        returns: List[AsyncResult] = []
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
            except (AttributeError, KeyError):  # pragma: no cover
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

        quma_result: str = self._access_quma(genomic_contents, read_contents)

        processed_quma: List[str] = self._process_raw_quma(quma_result)
        # Next, add this readname to the holding data frame
        int_df: pd.DataFrame = pd.DataFrame({position: processed_quma})
        return int_df

    @staticmethod
    def _access_quma(genomic_contents: str, read_contents: str) -> str:
        """Helper function to run internal QUMA tool.

        Args:

        Returns:
            str: output from quma command
        """

        result = quma.Quma(genomic_contents, read_contents).values

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
        """Config: pyllelic config object."""

        self.files_set: List[str] = files_set
        """List[str]: list of bam files analyzed."""

        with open(self.config.promoter_file, "r") as f:
            genome_base = f.readlines()
            genome_base_lines: List[str] = [s.rstrip("\n") for s in genome_base]
            self.genome_string: str = "".join(map(str, genome_base_lines))

        self._bam_output: Dict[str, BamOutput] = {}
        self._index_and_fetch()
        self.positions: List[str] = self._calculate_positions()
        """List[str]: list of genomic positions in the data."""

        self.cell_types = list(self._bam_output.keys())
        """List[str]: list of cell types in the data."""

        self.quma_results: Dict[str, QumaResult] = self._quma_full_threaded()
        """Dict[str, QumaResult]: list of QumaResults."""

        self.means: pd.DataFrame = self._process_means()
        """pd.DataFrame: dataframe of mean methylation values."""

        self.modes: pd.DataFrame = self._process_modes()
        """pd.DataFrame: dataframe of modes of methylation values."""

        self.diffs: pd.DataFrame = self._find_diffs(self.means, self.modes)
        """pd.DataFrame: dataframe of difference mean minus mode methylation values."""

        self.individual_data: pd.DataFrame = self._return_individual_data()
        """pd.DataFrame: dataframe of individual methylation values."""

    def _index_and_fetch(self) -> None:
        """Wrapper to call processing of each sam file.

        Args:
            genome_string (str): genomic sequence to align against

        Returns:
            list[str]: list of genomic positions analyzed
        """

        sam_path: List[Path] = [
            self.config.base_directory / "test" / f for f in self.files_set
        ]

        for sams in tqdm(sam_path, desc="Process BAM Files"):
            self._bam_output[str(sams)] = BamOutput(
                sams, self.genome_string, self.config
            )

    def _calculate_positions(self) -> List[str]:
        """Return sorted list of all positions analyzed.

        Returns:
            List[str]: sorted list of all positions analyzed.
        """

        positions = []
        for sublist in self._bam_output.values():
            for each in sublist.values.keys():
                positions.append(each)

        return sorted(list(set(positions)))

    def _quma_full_threaded(self) -> Dict[str, QumaResult]:
        """Run external QUMA methylation analysis on all specified cell lines,
        using multiprocessing library.

        Returns:
            Dict[str, QumaResult]: dictionary of cell_types and quma results
        """

        quma_results: Dict[str, QumaResult] = {}
        for name, bam_result in tqdm(
            self._bam_output.items(), desc="Process methylation"
        ):
            cell_line_name: str = name.split("_")[1]

            read_files: List[str] = [each for each in bam_result.values.values()]
            genomic_files: List[str] = [
                each for each in bam_result.genome_values.values()
            ]
            positions = bam_result.positions

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
            self.config.base_directory.joinpath(filename)
        )

        for name, each in self.quma_results.items():
            each.values.to_excel(writer, sheet_name=name)

        writer.save()

    def save_pickle(self, filename: str) -> None:
        """Save GenomicPositionData object as a pickled file."""

        with open(filename, "wb") as output_file:
            pickle.dump(self, output_file)

    @staticmethod
    def from_pickle(filename: str) -> GPD:
        """Read pickled GenomicPositionData back to an object."""

        with open(filename, "rb") as input_file:
            data = pickle.load(input_file)  # nosec

        return data

    def _process_means(self) -> pd.DataFrame:
        """Process the mean values at each position for each cell line.

        Returns:
            pd.DataFrame: dataframes of mean values for each position in each cell line
        """

        # Gives the means of each positions-- NOT the mean of the entire dataframe!
        working_df: pd.DataFrame = pd.DataFrame()
        for pos in self.positions:
            working_df[pos] = ""
            for key in self.quma_results.keys():
                values_list: List[float] = self._return_read_values(pos, key)
                if values_list:
                    pos_means: float = float(np.mean(values_list))
                else:  # No data or data doesn't meet minimums for analysis
                    pos_means = np.nan

                working_df.at[key, pos] = pos_means

        means_df: pd.DataFrame = working_df

        return means_df

    def _process_modes(self) -> pd.DataFrame:
        """Process the mode values at each position for each cell line.

        Returns:
            pd.DataFrame: dataframes of mode values for each position in each cell line
        """

        # Gives the modes of each positions-- NOT the mode of the entire dataframe!
        working_df: pd.DataFrame = pd.DataFrame()
        for pos in self.positions:
            working_df[pos] = ""
            for key in self.quma_results.keys():
                values_list: List[float] = self._return_read_values(pos, key)

                if values_list:
                    pos_modes: float = stats.mode(values_list)[0][0]
                else:  # No data or data doesn't meet minimums for analysis
                    pos_modes = np.nan

                working_df.at[key, pos] = pos_modes

        modes_df: pd.DataFrame = working_df

        return modes_df

    def _return_individual_data(self) -> pd.DataFrame:
        """Return a dataframe for methylation values at each position for each cell line.

        Returns:
            pd.DataFrame: methylation values for each position in each cell line
        """

        working_df: pd.DataFrame = pd.DataFrame()
        for pos in self.positions:
            working_df[pos] = ""  # Create position column in dataframe
            for key in self.quma_results.keys():
                values_list: List[float] = self._return_read_values(pos, key)
                if values_list:
                    data_for_df: Union[List[float], float] = values_list
                else:  # No data or data doesn't meet minimums for analysis
                    data_for_df = np.nan

                working_df.at[key, pos] = data_for_df

        return working_df

    def _return_read_values(
        self,
        pos: str,
        key: str,
        min_sites: int = MIN_NUMBER_READS,
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
            min_sites  # Only include is read has minimum methylation sites
        )

        values_list: List[float] = []
        # Check if this cell line has that position
        df: pd.DataFrame = self.quma_results[key].values
        if pos in df.columns:  # type: ignore[union-attr]

            values_to_check: List[str] = self._get_str_values(df, pos)
            for value in values_to_check:
                if not any(substring in value for substring in bad_values):
                    if len(value) >= min_num_meth_sites:
                        fraction_val: float = float(value.count("1")) / float(
                            len(value)
                        )  # number of methylated sites in each read
                        values_list.append(fraction_val)

        return values_list

    @staticmethod
    def _get_str_values(df: pd.DataFrame, pos: str) -> List[str]:
        """Return list of values a a position in given dataframe of methylation data.

        Args:
            df (pd.DataFrame): cell line dataframe
            pos (str): position being analyzed

        Returns:
            List[str]: list of values
        """
        return df.loc[:, pos].dropna().astype(str)  # type: ignore[union-attr]

    @staticmethod
    def _find_diffs(means_df: pd.DataFrame, modes_df: pd.DataFrame) -> pd.DataFrame:
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
    def _truncate_diffs(diffs_df: pd.DataFrame) -> pd.DataFrame:
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
            self.config.base_directory / str(filename + "_means.xlsx"),
            sheet_name="Means",
            index=True,
        )
        self.modes.to_excel(
            self.config.base_directory / str(filename + "_modes.xlsx"),
            sheet_name="Modes",
            index=True,
        )
        self.diffs.to_excel(
            self.config.base_directory / str(filename + "_diff.xlsx"),
            sheet_name="Diff",
            index=True,
        )

    @staticmethod
    def _create_histogram(
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

        fig: go.Figure = self._create_histogram(data, cell_line, position)
        fig.show()

    @staticmethod
    def _create_heatmap(
        df: pd.DataFrame, min_values: int, width: int, height: int, title_type: str
    ) -> go.Figure:
        """Generate a graph figure showing heatmap of mean methylation across
        cell lines.

        Args:
            df (pd.DataFrame): dataframe of mean methylation
            min_values (int): minimum number of points data must exist at a position
            width (int): figure width
            height (int): figure height
            title_type (str): type of figure being plotted

        Returns:
            go.Figure: plotly figure object
        """

        values_needed = len(df.index) - min_values
        df = df.loc[:, (df.isnull().sum(axis=0) <= values_needed)]

        fig: go.Figure = go.Figure()
        fig.add_trace(
            go.Heatmap(
                z=df,
                x=df.columns,
                y=df.index,
                hoverongaps=False,
            ),
        )
        fig.update_traces(
            dict(showscale=False, coloraxis=None, colorscale="RdBu_r"),
            selector={"type": "heatmap"},
        )

        fig.update_layout(
            title_text=f"{title_type} Methylation Heatmap",
            xaxis_title_text="Position",
            yaxis_title_text="Cell Line",
            # aspect="equal",
            template="seaborn",
            autosize=False,
            width=width,
            height=height,
        )

        return fig

    def heatmap(
        self,
        min_values: int,
        width: int = 800,
        height: int = 2000,
        cell_lines: Optional[List[str]] = None,
        data_type: str = "means",
    ) -> None:
        """Display a graph figure showing heatmap of mean methylation across
        cell lines.

        Args:
            min_values (int): minimum number of points data must exist at a position
            width (int): figure width, defaults to 800
            height (int): figure height, defaults to 2000
            cell_lines (Optional[List[str]]): set of cell lines to analyze,
            defaults to all cell lines.
            data_type (str): type of data to plot. Can to 'means', 'modes', or 'diffs'
        """

        title_type: str
        data: pd.DataFrame
        if data_type == "means":
            data = self.means
            title_type = "Mean"
        elif data_type == "modes":
            data = self.modes
            title_type = "Mode"
        elif data_type == "diffs":
            data = self.diffs
            title_type = "Difference"
        else:
            raise ValueError("Invalid data type")

        if cell_lines:
            data = data[data.index.isin(cell_lines)]

        fig: go.Figure = self._create_heatmap(
            data, min_values, width, height, title_type
        )
        fig.show()

    @staticmethod
    def _create_methylation_diffs_bar_graph(df: pd.DataFrame) -> go.Figure:
        """Generate a graph figure showing bar graph of significant methylation across
        cell lines.

        Args:
            df (pd.DataFrame): dataframe of significant methylation positions

        Returns:
            go.Figure: plotly figure object
        """
        data = df.pivot(index="position", columns="cellLine", values="ad_stat")
        data = data.dropna(axis=1, how="all").count(axis=1).to_frame()

        fig = go.Figure()
        fig.add_trace(go.Bar(x=data.index, y=data.values.flatten()))
        fig.update_layout(
            xaxis_type="linear",
            showlegend=False,
            title="Significant Methylation Differences",
            template="seaborn",
        )
        fig.update_xaxes(
            tickformat="r",
            tickangle=45,
            nticks=40,
            title="Position",
            range=[int(data.index.min()), int(data.index.max())],
        )
        fig.update_yaxes(title="# of significant differences")
        fig.update_traces(width=50)

        return fig

    def sig_methylation_differences(
        self, cell_lines: Optional[List[str]] = None
    ) -> None:
        """Display a graph figure showing a bar chart of significantly different
        mean / mode methylation across all or a subset of cell lines.

        Args:
            cell_lines (Optional[List[str]]): set of cell lines to analyze,
            defaults to all cell lines.
        """
        data = self.summarize_allelic_data(cell_lines)

        fig: go.Figure = self._create_methylation_diffs_bar_graph(data)
        fig.show()

    def generate_ad_stats(self) -> pd.DataFrame:
        """Generate Anderson-Darling normality statistics for an individual data df.

        Returns:
            pd.DataFrame: df of a-d test statistics
        """
        np.seterr(divide="ignore", invalid="ignore")  # ignore divide-by-zero errors
        df = self.individual_data
        df2 = df.applymap(self._anderson_darling_test)
        return df2

    def summarize_allelic_data(
        self, cell_lines: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """Create a dataframe only of likely allelic methylation positions.

        Args:
            cell_lines (Optional[List[str]]): set of cell lines to analyze,
            defaults to all cell lines.

        Returns:
            pd.DataFrame: dataframe of cell lines with likely allelic positions
        """

        if cell_lines:
            data = self.individual_data[self.individual_data.index.isin(cell_lines)]
        else:
            data = self.individual_data

        np.seterr(divide="ignore", invalid="ignore")  # ignore divide-by-zero errors
        sig_dict: Dict[str, List[np.ndarray]] = {
            "cellLine": [],
            "position": [],
            "ad_stat": [],
            "p_crit": [],
            "diff": [],
            "raw": [],
        }
        for index, row in data.iterrows():
            for column in row.index:
                value = row[column]
                if np.all(pd.notnull(value)):
                    good, stat, crits = self._anderson_darling_test(value)
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
    def _anderson_darling_test(raw_list: Optional[pd.Series]) -> AD_stats:
        """Run the Anderson-Darling normality test on methylation data at a point.

        Args:
            raw_list (pd.Series): list of fractional methylation levels per read.

        Returns:
            AD_stats[bool, float, List[np.ndarray]]: is the data significantly
                allelic (bool), A-D statistic (float), critical values (list)
        """

        if np.all(pd.notnull(raw_list)):
            stat: float
            crits: List[np.ndarray]
            stat, crits, _ = stats.anderson(raw_list)
            is_sig: bool = bool(stat > crits[4])
            return AD_stats(is_sig, stat, crits)
        return AD_stats(False, np.nan, [np.nan])


##################################################################################
##################################################################################

# Module level helper functions


def configure(
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

    config = Config(
        base_directory=Path(base_path),
        promoter_file=Path(base_path) / prom_file,
        results_directory=Path(base_path) / "results",
        bam_directory=Path(base_path) / "bam_output",
        analysis_directory=Path(base_path) / "test",
        promoter_start=prom_start,
        promoter_end=prom_end,
        chromosome=chrom,
        offset=offset,
    )

    return config


def make_list_of_bam_files(config: Config) -> List[str]:
    """Check analysis directory for all valid .bam files.

    Args:
        Config: pyllelic configuration options.

    Returns:
        list[str]: list of files
    """
    return [f.name for f in config.analysis_directory.iterdir() if f.suffix == ".bam"]
