#!/usr/bin/env python3
"""pytest unit tests for pyllelic."""

# Testing
import base64
import re
import unittest.mock as mock
from pathlib import Path
from typing import Any, Tuple, Union

import numpy as np

# Required libraries for test data
import pandas as pd
import pytest
from _pytest.monkeypatch import MonkeyPatch
from pytest_mock.plugin import MockerFixture

# Module to test
from pyllelic import pyllelic
from pyllelic.config import Config

from .inputs import (
    EXPECTED_ALLELIC_DATA,
    EXPECTED_BAM_INDIVIDUAL_POSITIONS,
    EXPECTED_BAM_OUTPUT_GENOME_VALUES,
    EXPECTED_BAM_OUTPUT_POSITIONS,
    EXPECTED_BAM_OUTPUT_VALUES,
    EXPECTED_INDIVIDUAL_POSITIONS,
    EXPECTED_INTERMEDIATE_DIFFS,
    EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA,
    EXPECTED_MEANS,
    EXPECTED_QUMA_VALUES,
    EXPECTED_STACKED_BAM,
    EXPECTED_WRITE_DF_OUTPUT,
    INPUT_READ_FILE,
    SAMPLE_BAI,
    SAMPLE_BAM,
    TEST_PROM_FILE,
)

# import os

PositionDataTuple = Tuple[Path, pyllelic.GenomicPositionData]


# Helper methods
@pytest.fixture(scope="session")
def set_up_genomic_position_data(
    tmp_path_factory: pytest.TempPathFactory,
) -> PositionDataTuple:
    tmp_path = tmp_path_factory.mktemp("data")
    p, _ = setup_bam_files(tmp_path)
    config = setup_config(p)
    INPUT_BAM_LIST = ["fh_test_tissue.bam"]
    return (
        p,
        pyllelic.GenomicPositionData(config=config, files_set=INPUT_BAM_LIST),
    )


@pytest.fixture(autouse=True)
def _mock_pool_apply_async(monkeypatch: MonkeyPatch) -> None:
    monkeypatch.setattr("multiprocessing.pool.Pool.apply_async", _mock_apply_async)


class MockPoolApplyResult:
    def __init__(self, func: Any, args: Any) -> None:
        self._func = func
        self._args = args

    def get(self, timeout: int = 0) -> Any:
        return self._func(*self._args)


def _mock_apply_async(
    self: Any,
    func: Any,
    args: Any = (),
    kwds: Any = None,
    callback: Any = None,
    error_callback: Any = None,
) -> MockPoolApplyResult:
    return MockPoolApplyResult(func, args)


def setup_bam_files(tmp_path: Path) -> Tuple[Path, Path]:
    d = tmp_path / "test"
    d.mkdir()
    fn_bam = "fh_test_tissue.bam"
    filepath_bam = d / fn_bam
    filepath_bam.write_bytes(base64.decodebytes(SAMPLE_BAM))
    fn_bai = "fh_test_tissue.bam.bai"
    filepath_bai = d / fn_bai
    filepath_bai.write_bytes(base64.decodebytes(SAMPLE_BAI))
    return tmp_path, filepath_bam


def setup_config(
    my_path: Path, offset: Union[int, None] = None, viz_backend: Union[str, None] = None
) -> Config:
    d = my_path
    prom_file = d / "test.txt"
    prom_file.write_text(TEST_PROM_FILE)
    TEST_START = 1293200
    TEST_END = 1296000
    TEST_CHR = "5"
    TEST_OFFSET = 1293000
    TEST_BACKEND = "plotly"
    if not offset:
        offset = TEST_OFFSET
    if not viz_backend:
        viz_backend = TEST_BACKEND
    config = pyllelic.configure(
        base_path=str(my_path),
        prom_file=str(prom_file),
        prom_start=TEST_START,
        prom_end=TEST_END,
        chrom=TEST_CHR,
        offset=offset,
        viz_backend=viz_backend,
    )

    return config


# Repeated test parameters
TEST_QUMA_RESULT = (
    "genome\t0\tATCGTAGTCGA\t2\t2,8\n"
    + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
    + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
    + "2\tquery2\tATCGATAGCATT\tATCGATAGCATT\tATCG-TAGTCGA\t"
    + "12\t5\t58.3\t1\t1\t0\t1\t100.0\t1A\t1\t1\n"
)


# Tests
def test_configure() -> None:
    """Test setting environment variables with mock object."""
    TEST_BASE_PATH = Path().cwd()
    TEST_PROM_FILE = TEST_BASE_PATH / "test.txt"
    TEST_START = "1"
    TEST_END = "2"
    TEST_CHR = "5"
    TEST_OFFSET = 0
    EXPECTED_RESULTS = TEST_BASE_PATH / "results"
    config = pyllelic.configure(
        base_path=str(TEST_BASE_PATH),
        prom_file=str(TEST_PROM_FILE),
        prom_start=int(TEST_START),
        prom_end=int(TEST_END),
        chrom=TEST_CHR,
        offset=TEST_OFFSET,
    )

    assert TEST_BASE_PATH == config.base_directory
    assert TEST_PROM_FILE == config.promoter_file
    assert int(TEST_START) == config.promoter_start
    assert int(TEST_END) == config.promoter_end
    assert EXPECTED_RESULTS == config.results_directory


def test_configure_fname_pattern() -> None:
    """Test setting environment variables with mock object."""
    TEST_BASE_PATH = Path().cwd()
    TEST_PROM_FILE = TEST_BASE_PATH / "test.txt"
    TEST_START = "1"
    TEST_END = "2"
    TEST_CHR = "5"
    TEST_OFFSET = 0
    EXPECTED_RESULTS = TEST_BASE_PATH / "results"
    TEST_FNAME_PATTERN = r"^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$"
    config = pyllelic.configure(
        base_path=str(TEST_BASE_PATH),
        prom_file=str(TEST_PROM_FILE),
        prom_start=int(TEST_START),
        prom_end=int(TEST_END),
        chrom=TEST_CHR,
        offset=TEST_OFFSET,
        fname_pattern=TEST_FNAME_PATTERN,
    )

    assert TEST_BASE_PATH == config.base_directory
    assert TEST_PROM_FILE == config.promoter_file
    assert int(TEST_START) == config.promoter_start
    assert int(TEST_END) == config.promoter_end
    assert EXPECTED_RESULTS == config.results_directory
    assert re.compile(TEST_FNAME_PATTERN) == config.fname_pattern


def test_make_list_of_bam_files(tmp_path: Path) -> None:
    config = setup_config(tmp_path)
    TEST_LIST = [
        Path("bad1.txt"),
        Path("bad2.bai"),
        Path("fh_good1_tissue.bam"),
        Path("fh_good2_tissue.bam"),
    ]
    EXPECTED = ["fh_good1_tissue.bam", "fh_good2_tissue.bam"]
    with mock.patch.object(
        pyllelic.Path,  # type:ignore[attr-defined]
        "iterdir",
    ) as mock_iterdir:
        mock_iterdir.return_value = TEST_LIST
        actual = pyllelic.make_list_of_bam_files(config)

    assert EXPECTED == actual


def test_pyllelic(tmp_path_factory: pytest.TempPathFactory) -> None:
    tmp_path = tmp_path_factory.mktemp("data")
    p, _ = setup_bam_files(tmp_path)
    config = setup_config(p)
    INPUT_BAM_LIST = ["fh_test_tissue.bam"]
    genomic_position_data = pyllelic.pyllelic(config=config, files_set=INPUT_BAM_LIST)
    print(repr(genomic_position_data.positions))

    assert genomic_position_data.positions == EXPECTED_BAM_INDIVIDUAL_POSITIONS
    assert genomic_position_data.cell_types == ["test"]


def test_compare_lines(
    set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
) -> None:
    _, genomic_position_data = set_up_genomic_position_data

    actual = pyllelic.compare_lines(genomic_position_data, "test", None, "means")

    assert actual.to_dict() == {"RMSE": {}}


# Tests of main classes


class Test_BamOutput:
    """Class to test BamOutput object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_bam_output(self, tmp_path: Path) -> pyllelic.BamOutput:
        p, fp_bam = setup_bam_files(tmp_path)
        config = setup_config(p)
        return pyllelic.BamOutput(
            sam_directory=fp_bam,
            genome_string=TEST_PROM_FILE,
            config=config,
        )

    def test_init(self, tmp_path: Path, set_up_bam_output: pyllelic.BamOutput) -> None:
        bam_output = set_up_bam_output

        assert bam_output.name == str(tmp_path / "test" / "fh_test_tissue.bam")
        assert bam_output.values == EXPECTED_BAM_OUTPUT_VALUES
        assert bam_output.positions == EXPECTED_BAM_OUTPUT_POSITIONS
        assert bam_output.genome_values == EXPECTED_BAM_OUTPUT_GENOME_VALUES

    def test_run_sam_and_extract_df(
        self, set_up_bam_output: pyllelic.BamOutput
    ) -> None:
        bam_output = set_up_bam_output
        actual_positions = bam_output._run_sam_and_extract_df(Path(bam_output.name))
        assert actual_positions == EXPECTED_BAM_OUTPUT_POSITIONS

    def test_write_bam_output(self, set_up_bam_output: pyllelic.BamOutput) -> None:
        bam_output = set_up_bam_output
        bam_output._write_bam_output(bam_output.positions, EXPECTED_STACKED_BAM)
        assert bam_output.values == EXPECTED_WRITE_DF_OUTPUT

    def test__genome_range(self, set_up_bam_output: pyllelic.BamOutput) -> None:
        bam_output = set_up_bam_output
        gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        result = bam_output._genome_range(position="8", genome_string=gen_str, offset=2)
        assert result == gen_str[6:36]
        assert isinstance(result, str)

    def test__genome_range_no_offset(
        self, set_up_bam_output: pyllelic.BamOutput
    ) -> None:
        bam_output = set_up_bam_output
        gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        new_config = setup_config(bam_output._config.base_directory, 2)
        bam_output._config = new_config
        result = bam_output._genome_range(position="8", genome_string=gen_str)
        assert result == gen_str[6:36]
        assert isinstance(result, str)

    def test__genome_range_invalid(self, set_up_bam_output: pyllelic.BamOutput) -> None:
        bam_output = set_up_bam_output
        gen_str = ""
        with pytest.raises(ValueError):
            _ = bam_output._genome_range(position="2", genome_string=gen_str, offset=40)

    def test_genome_parsing(self, set_up_bam_output: pyllelic.BamOutput) -> None:
        pass


class Test_QumaOutput:
    """Class to test QumaOutput object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_quma_output(self) -> pyllelic.QumaResult:
        INPUT_GENOMIC_FILE = ["CGGCGTAGGTAGGTTCGTACGAAGTCGTA"]

        return pyllelic.QumaResult(INPUT_READ_FILE, INPUT_GENOMIC_FILE, ["1293588"])

    def test_init(self, set_up_quma_output: pyllelic.QumaResult) -> None:
        quma_output = set_up_quma_output

        EXPECTED_QUMA_VALUES = pd.DataFrame.from_dict(
            {
                "1293588": {
                    0: "111",
                    1: "111",
                    2: "111",
                    3: "111",
                    4: "111",
                    5: "1111",
                }
            }
        )

        pd.testing.assert_frame_equal(quma_output.values, EXPECTED_QUMA_VALUES)

    def test_process_raw_quma(self, set_up_quma_output: pyllelic.QumaResult) -> None:
        quma_output = set_up_quma_output

        quma_reference = quma_output.quma_output[0].data

        EXPECTED = ["111", "111", "111", "111", "111", "1111"]
        actual = quma_output._process_raw_quma(quma_reference)
        assert EXPECTED == actual

    def test_process_raw_quma_below_min_alignment(
        self, set_up_quma_output: pyllelic.QumaResult
    ) -> None:
        quma_output = set_up_quma_output

        quma_reference = quma_output.quma_output[0].data
        quma_reference[0].res.perc = 10

        EXPECTED = ["FAIL", "111", "111", "111", "111", "1111"]
        actual = quma_output._process_raw_quma(quma_reference)
        assert EXPECTED == actual

    def test__pool_processing(self, set_up_quma_output: pyllelic.QumaResult) -> None:
        quma_output = set_up_quma_output

        EXPECTED_RESULTS = pd.DataFrame.from_dict(
            {
                "1293588": {
                    0: "111",
                    1: "111",
                    2: "111",
                    3: "111",
                    4: "111",
                    5: "1111",
                }
            }
        )

        actual = quma_output._pool_processing()

        pd.testing.assert_frame_equal(actual, EXPECTED_RESULTS)

    @mock.patch("pyllelic.pyllelic.signal.signal")
    def test__init_worker(
        self, mock_signal: mock.Mock, set_up_quma_output: pyllelic.QumaResult
    ) -> None:
        quma_output = set_up_quma_output
        quma_output._init_worker()
        mock_signal.assert_called_once()

    def test__thread_worker(self, set_up_quma_output: pyllelic.QumaResult) -> None:
        quma_output = set_up_quma_output
        TEST_READ_NAME = "1295094"
        TEST_GSEQ = ">genome\nATCGTAGTCGA"
        TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"

        EXPECTED: pd.DataFrame = pd.DataFrame({"1295094": ["1", "11"]})

        actual = quma_output._thread_worker(TEST_GSEQ, TEST_QSEQ, TEST_READ_NAME)

        pd.testing.assert_frame_equal(actual[0], EXPECTED)
        assert actual[1].values == EXPECTED_QUMA_VALUES

    def test_access_quma(self, set_up_quma_output: pyllelic.QumaResult) -> None:
        quma_output = set_up_quma_output
        TEST_GSEQ = ">genome\nATCGTAGTCGA"
        TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"

        actual = quma_output._access_quma(TEST_GSEQ, TEST_QSEQ)
        assert actual.values == EXPECTED_QUMA_VALUES


class Test_GenomicPositionData:
    """Class to test QumaOutput object initialization and read."""

    # pylint: disable=no-self-use

    def test_init(self, set_up_genomic_position_data: PositionDataTuple) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        assert genomic_position_data.positions == EXPECTED_BAM_INDIVIDUAL_POSITIONS
        assert genomic_position_data.cell_types == ["test"]

    # def test_save(self, set_up_genomic_position_data, mocker):
    #     _, genomic_position_data = set_up_genomic_position_data
    #     mocker.patch("pyllelic.pyllelic.pd")
    #     genomic_position_data.save()

    #     mocker.assert_called()

    def test_process_means(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        intermediate = EXPECTED_MEANS
        expected = intermediate.astype("object")

        np.array_equal(genomic_position_data.means.values, expected.values)
        # pd.testing.assert_frame_equal(genomic_position_data.means, expected)

    def test_histogram(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")
        TEST_POSITION = genomic_position_data.positions[0]
        TEST_CELL_LINE = genomic_position_data.means.index[0]

        genomic_position_data.histogram(TEST_CELL_LINE, TEST_POSITION)
        genomic_position_data.histogram(TEST_CELL_LINE, TEST_POSITION, backend="plotly")

        mocked_go.Figure.assert_called()
        mocked_go.Histogram.assert_called()

    def test_histogram_mpl(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_mpl = mocker.patch("pyllelic.visualization.sns")
        TEST_POSITION = genomic_position_data.positions[0]
        TEST_CELL_LINE = genomic_position_data.means.index[0]

        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="matplotlib"
        )
        genomic_position_data.config = new_config
        genomic_position_data.histogram(TEST_CELL_LINE, TEST_POSITION)

        mocked_mpl.histplot.assert_called_once()
        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="plotly"
        )
        genomic_position_data.config = new_config

    def test_histogram_error(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        TEST_POSITION = genomic_position_data.positions[0]
        TEST_CELL_LINE = genomic_position_data.means.index[0]

        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="other"
        )
        genomic_position_data.config = new_config
        with pytest.raises(ValueError):
            genomic_position_data.histogram(TEST_CELL_LINE, TEST_POSITION)
        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="plotly"
        )
        genomic_position_data.config = new_config

    def test_heatmap(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")

        genomic_position_data.heatmap(min_values=1)
        genomic_position_data.heatmap(min_values=1, backend="plotly")

        mocked_go.Figure.assert_called()
        mocked_go.Heatmap.assert_called()

    def test_heatmap_mpl(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_mpl = mocker.patch("pyllelic.visualization.sns")

        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="matplotlib"
        )
        genomic_position_data.config = new_config
        genomic_position_data.heatmap(min_values=1)

        mocked_mpl.heatmap.assert_called_once()
        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="plotly"
        )
        genomic_position_data.config = new_config

    def test_heatmap_error(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="other"
        )
        genomic_position_data.config = new_config
        with pytest.raises(ValueError):
            genomic_position_data.heatmap(min_values=1)
        new_config = setup_config(
            genomic_position_data.config.base_directory, viz_backend="plotly"
        )
        genomic_position_data.config = new_config

    def test_heatmap_cell_lines(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")

        TEST_CELL_LINE = "test.bam"
        genomic_position_data.heatmap(min_values=1, cell_lines=[TEST_CELL_LINE])

        mocked_go.Figure.assert_called_once()
        mocked_go.Heatmap.assert_called_once()

    def test_heatmap_modes(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")

        genomic_position_data.heatmap(min_values=1, data_type="modes")

        mocked_go.Figure.assert_called_once()
        mocked_go.Heatmap.assert_called_once()

    def test_heatmap_diffs(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")

        genomic_position_data.heatmap(min_values=1, data_type="diffs")

        mocked_go.Figure.assert_called_once()
        mocked_go.Heatmap.assert_called_once()

    def test_heatmap_allelic(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")

        genomic_position_data.heatmap(min_values=1, data_type="pvalue")

        mocked_go.Figure.assert_called_once()
        mocked_go.Heatmap.assert_called_once()

    def test_heatmap_invalid(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        with pytest.raises(ValueError):
            genomic_position_data.heatmap(min_values=1, data_type="FAKE")

    def test_reads_graph(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_px = mocker.patch("pyllelic.visualization.px")
        mocked_sp = mocker.patch("pyllelic.visualization.sp")
        mocked_df_plot = mocker.patch("pyllelic.visualization.pd.DataFrame.plot")
        mocked_mpl = mocker.patch("pyllelic.visualization.plt")
        mocked_mpl.subplots.return_value = (mocker.MagicMock(), mocker.MagicMock())

        genomic_position_data.reads_graph()
        genomic_position_data.reads_graph(backend="plotly")
        genomic_position_data.reads_graph(backend="matplotlib")

        mocked_px.bar.assert_called()
        mocked_sp.make_subplots.assert_called()
        mocked_df_plot.assert_called()
        mocked_mpl.subplots.assert_called()

    def test_reads_graph_cell_lines(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_px = mocker.patch("pyllelic.visualization.px")
        mocked_sp = mocker.patch("pyllelic.visualization.sp")
        CELL_LINES = ["test"]

        genomic_position_data.reads_graph(cell_lines=CELL_LINES)

        mocked_px.bar.assert_called()
        mocked_sp.make_subplots.assert_called_once()

    def test_reads_graph_too_many_cell_lines(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        CELL_LINES = ["test"]

        with pytest.raises(ValueError):
            genomic_position_data.reads_graph(cell_lines=CELL_LINES, max_graphs=0)

    def test_reads_graph_invalid(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        with pytest.raises(ValueError):
            genomic_position_data.reads_graph(backend="FAKE")

    def test_summarize_allelelic_data(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        with np.testing.suppress_warnings() as sup:  # ignore degrees of freedom warning
            sup.filter(RuntimeWarning, "Degrees of freedom")
            sup.filter(module=np.ma.core)
            CELL_LINES = ["test"]
            actual1 = genomic_position_data.summarize_allelic_data(
                cell_lines=CELL_LINES
            )
            actual2 = genomic_position_data.summarize_allelic_data()

        print("Actual")
        print(actual1.to_dict())
        print("Expected")
        print(EXPECTED_ALLELIC_DATA)
        # EXPECTED_ALLELIC_DATA.index = pd.RangeIndex(0, 0, 1)
        pd.testing.assert_frame_equal(EXPECTED_ALLELIC_DATA, actual1, atol=1e-2)
        pd.testing.assert_frame_equal(EXPECTED_ALLELIC_DATA, actual2, atol=1e-2)

    def test_sig_methylation_differences(
        self, set_up_genomic_position_data: PositionDataTuple, mocker: MockerFixture
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.visualization.go")
        mocked_mpl = mocker.patch("pyllelic.visualization.pd.DataFrame.plot")

        with np.testing.suppress_warnings() as sup:  # ignore deg of freedom warning
            sup.filter(RuntimeWarning, "Degrees of freedom")
            sup.filter(module=np.ma.core)
            genomic_position_data.sig_methylation_differences()
            genomic_position_data.sig_methylation_differences(backend="plotly")
            genomic_position_data.sig_methylation_differences(backend="matplotlib")

            mocked_go.Figure.assert_called()
            mocked_go.Bar.assert_called()
            mocked_mpl.assert_called_once()

    def test_sig_methylation_differences_invalid(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        with np.testing.suppress_warnings() as sup:  # ignore degrees of freedom warning
            sup.filter(RuntimeWarning, "Degrees of freedom")
            sup.filter(module=np.ma.core)
            with pytest.raises(ValueError):
                genomic_position_data.sig_methylation_differences(backend="FAKE")

    def test_anderson_darling_test_with_values(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        TEST_LIST = pd.Series([1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)])
        actual = genomic_position_data._anderson_darling_test(TEST_LIST)
        EXPECTED = (
            False,
            0.9267475729317516,
            np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
        )
        assert EXPECTED[0] == actual[0]
        np.testing.assert_almost_equal(EXPECTED[1], actual[1])  # type:ignore
        assert (EXPECTED[2] == actual[2]).all()

    def test_anderson_darling_test_with_bad(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        TEST_LIST = [np.nan]
        actual = genomic_position_data._anderson_darling_test(TEST_LIST)
        EXPECTED = (False, np.nan, [np.nan])
        assert EXPECTED[0] == actual[0]
        np.testing.assert_equal(EXPECTED[1], actual[1])
        assert EXPECTED[2] == actual[2]

    def test_generate_ad_stats(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        TEST_DF = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
        EXPECTED = pd.DataFrame.from_dict(
            {
                "TEST1": [
                    (
                        False,
                        0.9267475729317516,
                        np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
                    ),
                    (
                        False,
                        0.9267475729317516,
                        np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
                    ),
                    (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
                ],
                "TEST2": [
                    (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
                    (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
                    (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
                ],
                "TEST3": [
                    (False, np.nan, [np.nan]),
                    (False, np.nan, [np.nan]),
                    (
                        False,
                        0.7995458667216608,
                        np.array([0.72, 0.82, 0.984, 1.148, 1.365]),
                    ),
                ],
            },
            orient="index",
            columns=["1", "2", "3"],
        )
        genomic_position_data.individual_data = TEST_DF
        actual = genomic_position_data.generate_ad_stats()
        EXPECTED_TEST2_1 = EXPECTED.loc["TEST2", "1"]
        np.testing.assert_equal(EXPECTED_TEST2_1, actual.loc["TEST2", "1"])

    def test_write_means_modes_diffs(
        self, mocker: MockerFixture, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        TEST_FILENAME = "output"
        mocker.patch.object(genomic_position_data.means, "to_excel")
        mocker.patch.object(genomic_position_data.modes, "to_excel")
        mocker.patch.object(genomic_position_data.diffs, "to_excel")

        genomic_position_data.write_means_modes_diffs(TEST_FILENAME)

        genomic_position_data.means.to_excel.assert_called()
        genomic_position_data.modes.to_excel.assert_called()
        genomic_position_data.diffs.to_excel.assert_called()

    def test_save(
        self, mocker: MockerFixture, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        TEST_FILENAME = "output"
        # mocked_excel_writer = mocker.patch("pyllelic.pyllelic.pd.ExcelWriter")
        mocked_dataframe = mocker.patch("pyllelic.pyllelic.pd.DataFrame.to_excel")
        genomic_position_data.save(TEST_FILENAME)

        mocked_dataframe.assert_called()
        # mocked_excel_writer.ExcelWriter.assert_called()

    def test_save_pickle(
        self, mocker: MockerFixture, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data
        mocked_pickle = mocker.patch("pyllelic.pyllelic.pickle")

        TEST_FILENAME = "example.pickle"
        genomic_position_data.save_pickle(TEST_FILENAME)

        mocked_pickle.dump.assert_called()

    def test_from_pickle(self, mocker: MockerFixture) -> None:
        mocked_pickle = mocker.patch("pyllelic.pyllelic.pickle")

        TEST_FILENAME = "example.pickle"
        _ = pyllelic.GenomicPositionData.from_pickle(TEST_FILENAME)
        mocked_pickle.load.assert_called()

    def test__truncate_diffs(self) -> None:
        TEST_DIFFS = EXPECTED_INTERMEDIATE_DIFFS.astype("object")
        EXPECTED_TRUNCATED_DIFF = TEST_DIFFS.dropna(how="all")
        actual = pyllelic.GenomicPositionData._truncate_diffs(TEST_DIFFS)
        pd.testing.assert_frame_equal(actual, EXPECTED_TRUNCATED_DIFF)

    def test_return_individual_positions(
        self, set_up_genomic_position_data: PositionDataTuple
    ) -> None:
        _, genomic_position_data = set_up_genomic_position_data

        actual = genomic_position_data._return_individual_positions("test")

        pd.testing.assert_frame_equal(actual, EXPECTED_INDIVIDUAL_POSITIONS)
