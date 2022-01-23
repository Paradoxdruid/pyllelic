#!/usr/bin/env python3
"""pyllelic: a tool for detection of allelic-specific variation
   in reduced representation bisulfate DNA sequencing.
"""

import pickle  # nosec
import re
import signal
import warnings
from multiprocessing import Pool, cpu_count
from multiprocessing.pool import AsyncResult
from pathlib import Path
from typing import Dict, List, Match, NamedTuple, Optional, Pattern, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pysam
from numpy.typing import ArrayLike
from scipy import stats
from tqdm.auto import tqdm

from . import quma
from . import visualization as viz
from .config import Config

# Initialized multiprocessing limits
NUM_THREADS = cpu_count() - 1

# SANITY ALIGNMENT LIMITS
MIN_NUMBER_READS = 1
MIN_PERCENT_ALIGNMENT = 50


class AD_stats(NamedTuple):
    """Helper class for NamedTuple results from anderson_darling_test"""

    sig: bool
    stat: ArrayLike
    crits: List[ArrayLike]


class BamOutput:
    """Storage container to process BAM sequencing files and store processed results."""

    def __init__(self, sam_directory: Path, genome_string: str, config: Config) -> None:
        self._config: Config = config
        self.name: str = str(sam_directory)
        """str: path to bam file analyzed."""

        self.values: Dict[str, str] = {}
        """Dict[str, str]: dictionary of reads at a given position"""

        self.positions: List[str] = self._run_sam_and_extract_df(sam_directory)
        """"pd.Index: index of genomic positions in the bam file."""

        self.genome_values: Dict[str, str] = {}
        """Dict[str, str]: dictionary of read files and contents."""

        self._genome_parsing(genome_string)

    def _run_sam_and_extract_df(self, sams: Path) -> List[str]:
        """Process samfiles, pulling out sequence and position data
        and writing to folders/files.

        Args:
            sams (Path): path to a samfile

        Returns:
            List[str]: list of unique positions in the samfile
        """

        # Grab the promoter region of interest
        samm: pysam.AlignmentFile = pysam.AlignmentFile(str(sams), "rb")
        itern = samm.fetch(
            self._config.chromosome,
            self._config.promoter_start,
            self._config.promoter_end,
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

        return df2.index.unique().astype(str).tolist()  # type: ignore

    def _write_bam_output(self, positions: pd.Index, df: pd.Series) -> None:
        """Extract alignments from sequencing reads and output text strings
        in bam values dictionary.

        Args:
            positions (pd.Index): list of unique positions
            df (pd.Series): series of sequencing reads
        """

        for each1 in positions:
            alignments: List[str] = df.loc[each1].tolist()

            read_file: List[str] = []
            for index, each in enumerate(alignments):
                read_file.append(str(">read" + str(index)))
                read_file.append(each)

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

        Raises:
            ValueError: incorrect position
        """

        if not offset:
            offset = self._config.offset

        start: int = int(position) - offset
        end: int = (int(position) - offset) + 30
        if genome_string[start:end]:
            return genome_string[start:end]
        raise ValueError("Position not in genomic promoter file")

    def _genome_parsing(self, genome_string: str) -> None:
        """Writes out a list of genomic sequence strings for comparison to read data.

        Args:
            genome_string (str): genomic sequence
        """

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
        self.quma_output: List[quma.Quma] = []
        """List[quma.Quma]: list of Quma result objects."""

        self.values: pd.DataFrame = self._pool_processing()
        """pd.DataFrame: dataframe of quma methylation analysis values."""

    @staticmethod
    def _process_raw_quma(quma_data: List[quma.Reference]) -> List[str]:
        """Spit and process raw quma results into a pandas dataframe.

        Args:
            quma_data (List[quma.Reference]): list of quma References.

        Returns:
            pd.DataFrame: intermediate dataframe to append
        """

        dots = [
            res.res.val if res.res.perc > MIN_PERCENT_ALIGNMENT else "FAIL"
            for res in quma_data
        ]

        return sorted(dots, key=lambda t: t.count("1"))

    def _pool_processing(self) -> pd.DataFrame:
        """Helper to run quma in a multiprocessing pool.

        Returns:
            pd.DataFrame: dataframe of quma results
        """

        # Set up a holding data frame from all the data
        holding_df: pd.DataFrame = pd.DataFrame()

        returns: List[AsyncResult[Tuple[pd.DataFrame, quma.Quma]]] = []
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
        for int_df, quma_obj in results:
            try:
                holding_df = pd.concat([holding_df, int_df], axis=1)
                self.quma_output.append(quma_obj)
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
    ) -> Tuple[pd.DataFrame, quma.Quma]:
        """Queue worker for quma functions.

        Args:
            genomic_contents (str): genomic sequence
            read_contents (str): query sequence
            position(str): position of reads

        Returns:
            Tuple[pd.DataFrame, quma.Quma]: dataframe of quma results and raw quma objs
        """

        quma_obj: quma.Quma = self._access_quma(genomic_contents, read_contents)

        processed_quma: List[str] = self._process_raw_quma(quma_obj.data)
        # Next, add this readname to the holding data frame
        int_df: pd.DataFrame = pd.DataFrame({position: processed_quma})
        return int_df, quma_obj

    @staticmethod
    def _access_quma(genomic_contents: str, read_contents: str) -> quma.Quma:
        """Helper function to run internal QUMA tool.

        Args:
            genomic_contents (str): genome sequence
            read_contents (str): query sequence

        Returns:
            quma.Quma: quma object
        """

        return quma.Quma(genomic_contents, read_contents)


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

        self._bam_output: Dict[str, BamOutput] = self._index_and_fetch()

        self.positions: List[str] = self._calculate_positions()
        """List[str]: list of genomic positions in the data."""

        self.file_names: List[str] = list(self._bam_output.keys())
        """List[str]: list of bam filenames in the data."""

        self.cell_types: List[str] = self._retrieve_cell_names()
        """List[str]: list of cell types in the data."""

        self.quma_results: Dict[str, QumaResult] = self._quma_full_threaded()
        """Dict[str, QumaResult]: list of QumaResults."""

        with tqdm(total=4, desc="Summary Statistics") as pbar:
            self.means: pd.DataFrame = self._process_means()
            """pd.DataFrame: dataframe of mean methylation values."""
            pbar.update(1)

            self.modes: pd.DataFrame = self._process_modes()
            """pd.DataFrame: dataframe of modes of methylation values."""
            pbar.update(1)

            self.diffs: pd.DataFrame = self._find_diffs(self.means, self.modes)
            """pd.DataFrame: df of difference mean minus mode methylation values."""
            pbar.update(1)

            self.individual_data: pd.DataFrame = self._return_individual_data()
            """pd.DataFrame: dataframe of individual methylation values."""

            self.positions = self.means.columns.tolist()
            pbar.update(1)

    def _index_and_fetch(self) -> Dict[str, BamOutput]:
        """Wrapper to call processing of each sam file.

        Returns:
            Dict[str, BamOutput]: dictionary of bam outputs"""

        sam_path: List[Path] = [
            self.config.base_directory / "test" / f for f in self.files_set
        ]
        bam_dict: Dict[str, BamOutput] = {}
        for sams in tqdm(sam_path, desc="Process BAM Files"):
            bam_dict[str(sams)] = BamOutput(sams, self.genome_string, self.config)
        return bam_dict

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

    def _retrieve_cell_names(self) -> List[str]:
        """Retrieve shortened cell line names from the list of filenames.

        Returns:
            List[str]: list of cell line names.
        """
        cells: List[str] = []
        for each in self.file_names:
            match: Optional[Match[str]] = re.match(
                self.config.fname_pattern, Path(each).name
            )
            cell_line_name: str
            if isinstance(match, Match):
                cell_line_name = match.group(1)
            else:  # pragma: no cover
                cell_line_name = Path(each).name
            cells.append(cell_line_name)
        return cells

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
            match: Optional[Match[str]] = re.match(
                self.config.fname_pattern, Path(name).name
            )
            cell_line_name: str
            if isinstance(match, Match):
                cell_line_name = match.group(1)
            else:  # pragma: no cover
                cell_line_name = Path(name).name

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
            filename (str): Filename to save to. Defaults to "output.xlsx".
        """

        writer: pd.ExcelWriter = pd.ExcelWriter(
            self.config.base_directory.joinpath(filename)
        )

        for name, each in self.quma_results.items():
            each.values.to_excel(writer, sheet_name=name)

        writer.save()

    def save_pickle(self, filename: str) -> None:
        """Save GenomicPositionData object as a pickled file.

        Args:
            filename (str): filename to save pickle
        """

        with open(filename, "wb") as output_file:
            pickle.dump(self, output_file)

    @staticmethod
    def from_pickle(filename: str) -> "GenomicPositionData":
        """Read pickled GenomicPositionData back to an object.

        Args:
            filename (str): filename to read pickle

        Returns:
            GenomicPositionData: GenomicPositionData object
        """

        with open(filename, "rb") as input_file:
            data: GenomicPositionData = pickle.load(input_file)  # nosec

        return data

    def _process_means(self) -> pd.DataFrame:
        """Process the mean values at each individual position for each cell line.

        Returns:
            pd.DataFrame: dataframes of mean values for each position in each cell line.
        """

        df_list = {
            cell: self.return_individual_positions(cell).apply(pd.to_numeric).mean()
            for cell in self.cell_types
        }

        df = pd.DataFrame.from_dict(df_list).T
        df.columns = df.columns.astype(str)
        return df

    def _process_modes(self) -> pd.DataFrame:
        """Process the mode at each individual position for each cell line.

        Returns:
            pd.DataFrame: dataframes of mode values for each individual position.
        """

        df_list = {
            cell: self.return_individual_positions(cell)
            .apply(pd.to_numeric)
            .mode()
            .iloc[0]
            for cell in self.cell_types
        }

        df = pd.DataFrame.from_dict(df_list).T
        df.columns = df.columns.astype(str)
        return df

    def _return_individual_data(self) -> pd.DataFrame:
        """Return a dataframe for methylation values at each position for each cell line.

        Returns:
            pd.DataFrame: methylation values for each position in each cell line
        """

        df_list: List[pd.DataFrame] = []
        for cell in self.cell_types:
            ind_pos: pd.DataFrame = self.return_individual_positions(cell).apply(
                pd.to_numeric
            )
            ind_pos = ind_pos.stack().reset_index(level=0, drop=True)
            ind_pos = (
                ind_pos.groupby(ind_pos.index)
                .apply(list)
                .to_frame()
                .transpose()
                .rename(index={0: cell})
            )

            df_list.append(ind_pos)

        df = pd.concat(df_list)
        df.columns = df.columns.astype(str)
        return df

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

    def histogram(
        self, cell_line: str, position: str, backend: Optional[str] = None
    ) -> None:
        """Display a graph figure showing fractional methylation in
        a given cell line at a given site.

        Args:
            cell_line (str): name of cell line
            position (str): genomic position
            backend (Optional[str]): plotting backend to override default

        Raises:
            ValueError: invalid plotting backend
            ValueError: No data available at that position
        """
        data = self.individual_data
        if not backend:
            backend = self.config.viz_backend

        if np.all(np.isnan(data.loc[cell_line, position])):
            raise ValueError("No data available at that position")

        if backend == "plotly":
            fig: go.Figure = viz._create_histogram(data, cell_line, position, backend)
            fig.show()
            return
        if backend == "matplotlib":
            _: plt.Figure = viz._create_histogram(data, cell_line, position, backend)
            plt.show()
            return
        raise ValueError("Invalid plotting backend")

    def reads_graph(
        self, cell_lines: Optional[List[str]] = None, backend: Optional[str] = None
    ) -> None:
        """Display a graph figure showing methylation of reads across cell lines.

        Args:
            cell_lines (Optional[List[str]]): set of cell lines to analyze,
                                            defaults to all cell lines.
            backend (Optional[str]): plotting backend to override default

        Raises:
            ValueError: invalid plotting backend
        """

        data = self.individual_data
        if not backend:
            backend = self.config.viz_backend

        if cell_lines:
            data = data[data.index.isin(cell_lines)]

        if backend == "plotly":
            fig: go.Figure = viz._make_stacked_fig(data, backend)
            fig.show()
            return
        if backend == "matplotlib":
            _: plt.Figure = viz._make_stacked_fig(data, backend)
            plt.show()
            return
        raise ValueError("Invalid plotting backend")

    def heatmap(
        self,
        min_values: int,
        width: int = 800,
        height: int = 2000,
        cell_lines: Optional[List[str]] = None,
        data_type: str = "means",
        backend: Optional[str] = None,
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
            backend (Optional[str]): plotting backend to override default

        Raises:
            ValueError: invalid data type
            ValueError: invalid plotting backend
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

        if not backend:
            backend = self.config.viz_backend

        if backend == "plotly":
            fig: go.Figure = viz._create_heatmap(
                data, min_values, width, height, title_type, backend
            )
            fig.show()
            return

        if backend == "matplotlib":
            _: plt.Figure = viz._create_heatmap(
                data, min_values, width, height, title_type, backend
            )
            plt.show()
            return

        raise ValueError("Invalid plotting backend")

    def sig_methylation_differences(
        self,
        cell_lines: Optional[List[str]] = None,
        backend: Optional[str] = None,
    ) -> None:
        """Display a graph figure showing a bar chart of significantly different
        mean / mode methylation across all or a subset of cell lines.

        Args:
            cell_lines (Optional[List[str]]): set of cell lines to analyze,
                                            defaults to all cell lines.
            backend (Optional[str]): plotting backend to override default

        Raises:
            ValueError: invalid plotting backend
        """
        data = self.summarize_allelic_data(cell_lines)

        if not backend:
            backend = self.config.viz_backend

        if backend == "plotly":
            fig: go.Figure = viz._create_methylation_diffs_bar_graph(data, backend)
            fig.show()
            return
        if backend == "matplotlib":
            _: plt.Figure = viz._create_methylation_diffs_bar_graph(data, backend)
            plt.show()
            return
        raise ValueError("Invalid plotting backend")

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
        sig_dict: Dict[str, List[ArrayLike]] = {
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
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
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
            AD_stats: is the data significantly
                allelic (bool), A-D statistic (float), critical values (list)
        """

        if np.all(pd.notnull(raw_list)):
            stat: ArrayLike
            crits: List[ArrayLike]
            stat, crits, _ = stats.anderson(raw_list)
            is_sig: bool = bool(stat > crits[4])  # type: ignore
            return AD_stats(is_sig, stat, crits)
        return AD_stats(False, np.nan, [np.nan])

    def return_individual_positions(self, cell_line: str) -> pd.DataFrame:
        """Return a dataframe of methylation at genomic position resolution.

        Args:
            cell_line (str): Cell line to analyze.

        Returns:
            pd.DataFrame: dataframe of methylation at genomic position resolution.
        """

        class PointData(NamedTuple):
            """Named tuple to cointain individual Quma point data."""

            cell: str
            read_position: int
            read: str
            positions: List[int]
            values: List[str]

        def _find_positions(
            cell: str, pos: str, data_list: List[quma.Reference], gseq: str
        ) -> List[PointData]:
            """Find individual resolution positions and data in quma results.

            Args:
                cell (str): cell line.
                pos (str): position of read.
                data_list (List[quma.Reference]): list of quma read data
                gseq (str): reference genomic sequence of read position

            Returns:
                List[PointData]: list of data points from quma results.
            """
            output: List[PointData] = []
            for each in data_list:
                read: str = each.fasta.com
                offset: List[int] = [
                    m.start() for m in re.finditer(each.res.gAli, gseq)
                ]
                real_off_set: int
                if offset:
                    real_off_set = offset[0]
                else:
                    real_off_set = 0
                init_pos_list: List[int] = [
                    m.start() for m in re.finditer("CG", each.res.gAli)
                ]
                pos_list: List[int] = [
                    (this + int(pos) + real_off_set) for this in init_pos_list
                ]

                my_point: PointData = PointData(
                    cell=cell,
                    read_position=int(pos),
                    read=read,
                    positions=pos_list,
                    values=list(each.res.val),
                )
                output.append(my_point)

            return output

        def _collate_positions(output: List[PointData]) -> Dict[int, List[str]]:
            """Collect together list of quma data points into position-indexed dict.

            Args:
                output (List[PointData]): list of quma read data points

            Returns:
                Dict[int, str]: position-indexed genomic resolution methylation dict.
            """
            result_dict: Dict[int, List[str]] = {}
            for each in output:
                for idx, val in enumerate(each.positions):
                    try:
                        my_val: str = each.values[idx]
                    except IndexError:  # pragma: no cover
                        pass
                    if my_val in ["1", "0"]:
                        if result_dict.get(val):
                            result_dict[val].append(my_val)
                        else:
                            result_dict[val] = [my_val]

            result_dict = {k: sorted(v) for k, v in result_dict.items()}
            return result_dict

        def _make_pos_df(data: Dict[int, List[str]]) -> pd.DataFrame:
            """Convert dictionary of genomic resolution data into dataframe.

            Args:
                data (Dict[int, List[str]]): pos-indexed high res. methylation dict.

            Returns:
                pd.DataFrame: dataframe of pos-indexed high resolution methylation.
            """

            int_df = pd.DataFrame.from_dict(data, orient="index").transpose()
            return int_df

        query: QumaResult = self.quma_results[cell_line]
        positions: List[str] = query._positions
        data: List[quma.Quma] = query.quma_output
        big_list: List[List[PointData]] = []
        for pos, dat in zip(positions, data):
            big_list.append(_find_positions(cell_line, pos, dat.data, dat._gseq))

        final_list: List[PointData] = [item for sublist in big_list for item in sublist]
        full_data: Dict[int, List[str]] = _collate_positions(final_list)

        final_df: pd.DataFrame = _make_pos_df(full_data)
        return final_df


##################################################################################
##################################################################################

# Module level helper functions


def configure(
    base_path: str,
    prom_file: str,
    prom_start: int,
    prom_end: int,
    chrom: str,
    offset: int,
    test_dir: Optional[str] = None,
    fname_pattern: Optional[str] = None,
    viz_backend: Optional[str] = None,
    results_dir: Optional[str] = None,
) -> Config:
    """Helper method to set up all our environmental variables, such as for testing.

    Args:
        base_path (str): directory where all processing will occur, put .bam files
                         in "test" sub-directory in this folder
        prom_file (str): filename of genmic sequence of promoter region of interest
        prom_start (int): start position to analyze in promoter region
        prom_end (int): final position to analyze in promoter region
        chrom (str): chromosome promoter is located on
        offset (int): genomic position of promoter to offset reads
        test_dir(Optional[str]): name of test directory where bam files are located
        fname_pattern(Optional[str]): regex pattern for processing filenames
        viz_backend (Optional[str]): which plotting backend to use
        results_dir (Optional[str]): name of results directory

    Returns:
        Config: configuration dataclass instance.
    """

    if not test_dir:
        test_dir = "test"

    fname_pattern_compiled: Pattern[str]
    if not fname_pattern:
        fname_pattern_compiled = re.compile(r"^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$")
    else:
        fname_pattern_compiled = re.compile(fname_pattern)

    if not viz_backend:
        viz_backend = "plotly"

    if not results_dir:
        results_dir = "results"

    config = Config(
        base_directory=Path(base_path),
        promoter_file=Path(base_path) / prom_file,
        results_directory=Path(base_path) / results_dir,
        analysis_directory=Path(base_path) / test_dir,
        promoter_start=prom_start,
        promoter_end=prom_end,
        chromosome=chrom,
        offset=offset,
        viz_backend=viz_backend,
        fname_pattern=fname_pattern_compiled,
    )

    return config


def make_list_of_bam_files(config: Config) -> List[str]:
    """Check analysis directory for all valid .bam files.

    Args:
        config (Config): pyllelic configuration options.

    Returns:
        list[str]: list of files
    """
    return [f.name for f in config.analysis_directory.iterdir() if f.suffix == ".bam"]


def pyllelic(config: Config, files_set: List[str]) -> GenomicPositionData:
    """Wrapper to call pyllelic routines.

    Args:
        config (Config): pyllelic config object.
        files_set (List[str]): list of bam files to analyze.

    Returns:
        GenomicPositionData: GenomicPositionData pyllelic object.
    """
    return GenomicPositionData(config=config, files_set=files_set)
