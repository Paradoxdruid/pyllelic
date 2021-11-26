#!/usr/bin/env python3
"""pytest unit tests for pyllelic."""

# Testing
import base64
import pytest
import unittest.mock as mock
import os

from inputs import (
    SAMPLE_BAM,
    SAMPLE_BAI,
    TEST_PROM_FILE,
    EXPECTED_BAM_OUTPUT_POSITIONS,
    EXPECTED_BAM_OUTPUT_VALUES,
    EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA,
    EXPECTED_BAM_OUTPUT_GENOME_VALUES,
    EXPECTED_STACKED_BAM,
    EXPECTED_WRITE_DF_OUTPUT,
    INPUT_READ_FILE,
    EXPECTED_MEANS,
)

# Required libraries for test data
import pandas as pd

import numpy as np
from pathlib import Path

# Module to test
from pyllelic import pyllelic


# Helper methods
@pytest.fixture(scope="session")
def set_up_genomic_position_data(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("data")
    p, _ = setup_bam_files(tmp_path)
    config = setup_config(p)
    INPUT_BAM_LIST = ["fh_test.bam"]
    return (
        p,
        pyllelic.GenomicPositionData(config=config, files_set=INPUT_BAM_LIST),
    )


@pytest.fixture(autouse=True)
def mock_pool_apply_async(monkeypatch):
    monkeypatch.setattr("multiprocessing.pool.Pool.apply_async", _mock_apply_async)


def _mock_apply_async(
    self, func, args=(), kwds=None, callback=None, error_callback=None
):
    return MockPoolApplyResult(func, args)


class MockPoolApplyResult:
    def __init__(self, func, args):
        self._func = func
        self._args = args

    def get(self, timeout=0):
        return self._func(*self._args)


def setup_bam_files(tmp_path):
    d = tmp_path / "test"
    d.mkdir()
    fn_bam = "fh_test.bam"
    filepath_bam = d / fn_bam
    filepath_bam.write_bytes(base64.decodebytes(SAMPLE_BAM))
    fn_bai = "fh_test.bam.bai"
    filepath_bai = d / fn_bai
    filepath_bai.write_bytes(base64.decodebytes(SAMPLE_BAI))
    return tmp_path, filepath_bam


def setup_config(my_path):
    d = my_path
    prom_file = d / "test.txt"
    prom_file.write_text(TEST_PROM_FILE)
    TEST_START = "1293200"
    TEST_END = "1296000"
    TEST_CHR = "5"
    TEST_OFFSET = 1298163
    config = pyllelic.configure(
        base_path=my_path,
        prom_file=prom_file,
        prom_start=TEST_START,
        prom_end=TEST_END,
        chrom=TEST_CHR,
        offset=TEST_OFFSET,
    )

    return config


# Tests
def test_configure():
    """Test setting environment variables with mock object."""
    TEST_BASE_PATH = Path().cwd()
    TEST_PROM_FILE = TEST_BASE_PATH / "test.txt"
    TEST_START = "1"
    TEST_END = "2"
    TEST_CHR = "5"
    TEST_OFFSET = 0
    EXPECTED_RESULTS = TEST_BASE_PATH / "results"
    config = pyllelic.configure(
        base_path=TEST_BASE_PATH,
        prom_file=TEST_PROM_FILE,
        prom_start=TEST_START,
        prom_end=TEST_END,
        chrom=TEST_CHR,
        offset=TEST_OFFSET,
    )

    assert TEST_BASE_PATH == config.base_directory
    assert TEST_PROM_FILE == config.promoter_file
    assert TEST_START == config.promoter_start
    assert TEST_END == config.promoter_end
    assert EXPECTED_RESULTS == config.results_directory


def test_make_list_of_bam_files(tmp_path):
    """Test making list of bam files, mocking config."""
    config = setup_config(tmp_path)
    TEST_LIST = [
        Path("bad1.txt"),
        Path("bad2.bai"),
        Path("good1.bam"),
        Path("good2.bam"),
    ]
    EXPECTED = ["good1.bam", "good2.bam"]
    with mock.patch.object(pyllelic.Path, "iterdir") as mock_iterdir:
        mock_iterdir.return_value = TEST_LIST
        actual = pyllelic.make_list_of_bam_files(config)

    assert EXPECTED == actual


# Tests of main classes


class Test_BamOutput:
    """Class to test BamOutput object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_bam_output(self, tmp_path):
        p, fp_bam = setup_bam_files(tmp_path)
        config = setup_config(p)
        return pyllelic.BamOutput(
            sam_directory=fp_bam,
            genome_string=TEST_PROM_FILE,
            config=config,
        )

    def test_init(self, tmp_path, set_up_bam_output):
        bam_output = set_up_bam_output

        assert bam_output.name == str(tmp_path / "test" / "fh_test.bam")
        assert bam_output.values == EXPECTED_BAM_OUTPUT_VALUES
        assert bam_output.positions == EXPECTED_BAM_OUTPUT_POSITIONS
        assert bam_output.genome_values == EXPECTED_BAM_OUTPUT_GENOME_VALUES

    def test_run_sam_and_extract_df(self, set_up_bam_output):
        bam_output = set_up_bam_output
        actual_positions = bam_output._run_sam_and_extract_df(Path(bam_output.name))
        assert actual_positions == EXPECTED_BAM_OUTPUT_POSITIONS

    def test_write_bam_output(self, set_up_bam_output):
        bam_output = set_up_bam_output
        # print(f"Initial bam values:\n{bam_output.values}")
        bam_output._write_bam_output(bam_output.positions, EXPECTED_STACKED_BAM)
        # print(f"Revised bam values:\n{bam_output.values}")
        assert bam_output.values == EXPECTED_WRITE_DF_OUTPUT

    @mock.patch("pyllelic.pyllelic.pysam")
    def test_pysam_index(self, mock_pysam, set_up_bam_output):
        bam_output = set_up_bam_output
        TEST_PATH = Path().cwd()
        bam_output._pysam_index(TEST_PATH)
        mock_pysam.index.assert_called_once_with(os.fspath(TEST_PATH))

    def test__genome_range(self, set_up_bam_output):
        """Check if correct genome string is returned."""
        bam_output = set_up_bam_output
        gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        result = bam_output._genome_range(position=2, genome_string=gen_str, offset=40)
        assert result == gen_str[8:37]
        assert isinstance(result, str)

    def test__genome_range_no_offset(self, set_up_bam_output):
        """Check if correct genome string is returned with no offset given."""
        bam_output = set_up_bam_output
        gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        bam_output._config.offset = 40
        result = bam_output._genome_range(position=2, genome_string=gen_str)
        assert result == gen_str[8:37]
        assert isinstance(result, str)

    def test_genome_parsing(self, set_up_bam_output):
        pass


class Test_QumaOutput:
    """Class to test QumaOutput object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_quma_output(self):

        INPUT_GENOMIC_FILE = ["CGGCGTAGGTAGGTTCGTACGAAGTCGTA"]

        return pyllelic.QumaResult(INPUT_READ_FILE, INPUT_GENOMIC_FILE, ["1293588"])

    def test_init(self, set_up_quma_output):
        quma_output = set_up_quma_output

        EXPECTED_QUMA_VALUES = pd.DataFrame.from_dict(
            {
                "1293588": {
                    0: "G11-1",
                    1: "G11-1",
                    2: "G11-1",
                    3: "G11-1",
                    4: "G11-1",
                    5: "1111",
                }
            }
        )

        pd.testing.assert_frame_equal(quma_output.values, EXPECTED_QUMA_VALUES)

    def test_process_raw_quma(self, set_up_quma_output):
        quma_output = set_up_quma_output
        TEST_QUMA_RESULT = (
            "genome\t0\tATCGTAGTCGA\t2\t2,8\n"
            + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
            + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
            + "2\tquery2\tATCGATAGCATT\tATCGATAGCATT\tATCG-TAGTCGA\t"
            + "12\t5\t58.3\t1\t1\t0\t1\t100.0\t1A\t1\t1\n"
        )
        EXPECTED = ["1A", "11"]
        actual = quma_output._process_raw_quma(TEST_QUMA_RESULT)
        assert EXPECTED == actual

    def test__pool_processing(self, set_up_quma_output):
        quma_output = set_up_quma_output

        EXPECTED_RESULTS = pd.DataFrame.from_dict(
            {
                "1293588": {
                    0: "G11-1",
                    1: "G11-1",
                    2: "G11-1",
                    3: "G11-1",
                    4: "G11-1",
                    5: "1111",
                }
            }
        )

        actual = quma_output._pool_processing()

        print(actual.to_dict())
        pd.testing.assert_frame_equal(actual, EXPECTED_RESULTS)

    @mock.patch("pyllelic.pyllelic.signal.signal")
    def test__init_worker(self, mock_signal, set_up_quma_output):
        quma_output = set_up_quma_output
        quma_output._init_worker()
        mock_signal.assert_called_once()

    def test__thread_worker(self, set_up_quma_output):
        quma_output = set_up_quma_output
        TEST_READ_NAME = "1295094"
        TEST_GSEQ = ">genome\nATCGTAGTCGA"
        TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"

        EXPECTED: pd.DataFrame = pd.DataFrame({"1295094": ["1", "11"]})

        actual = quma_output._thread_worker(TEST_GSEQ, TEST_QSEQ, TEST_READ_NAME)

        pd.testing.assert_frame_equal(actual, EXPECTED)

    def test_access_quma(self, set_up_quma_output):
        quma_output = set_up_quma_output
        TEST_GSEQ = ">genome\nATCGTAGTCGA"
        TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"

        EXPECTED = (
            "genome\t0\tATCGTAGTCGA\t1\t0\n1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\t"
            + "ATCGTAGTCGA\t11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n2\tquery2\t"
            + "ATCGATAGCATT\tATCG-TAGT\tATCGATAGC\t9\t1\t88.9\t1\t1\t0\t1\t"
            + "100.0\t1\t1\t1\n"
        )
        actual = quma_output._access_quma(TEST_GSEQ, TEST_QSEQ)
        assert EXPECTED == actual


class Test_GenomicPositionData:
    """Class to test QumaOutput object initialization and read."""

    # pylint: disable=no-self-use

    def test_init(self, set_up_genomic_position_data):
        p, genomic_position_data = set_up_genomic_position_data
        positions = []
        for each in EXPECTED_BAM_OUTPUT_POSITIONS:
            positions.append(each)

        EXPECTED_POSITIONS = sorted(set(positions))

        assert genomic_position_data.positions == EXPECTED_POSITIONS
        assert genomic_position_data.cell_types == [str(p / "test" / "fh_test.bam")]

    def test_process_means(self, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data

        intermediate = EXPECTED_MEANS
        expected = intermediate.astype("object")

        pd.testing.assert_frame_equal(genomic_position_data.means, expected)

    def test__create_histogram(self, mocker, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.pyllelic.go")
        TEST_CELL_LINE = "TEST1"
        TEST_POSITION = "1"
        intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
        TEST_DATA = intermediate.astype("object")

        _ = genomic_position_data._create_histogram(
            TEST_DATA, TEST_CELL_LINE, TEST_POSITION
        )

        mocked_go.Figure.assert_called_once()
        mocked_go.Histogram.assert_called_once_with(
            x=TEST_DATA.loc[TEST_CELL_LINE, TEST_POSITION],
            xbins=dict(
                start=-0.1,
                end=1.1,
                size=0.2,
            ),
        )

    def test_histogram(self, set_up_genomic_position_data, mocker):
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.pyllelic.go")
        TEST_CELL_LINE = "test.bam"
        TEST_POSITION = genomic_position_data.positions[0]

        genomic_position_data.histogram(TEST_CELL_LINE, TEST_POSITION)

        mocked_go.Figure.assert_called_once()
        mocked_go.Histogram.assert_called_once()

    def test__create_heatmap(self, mocker, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.pyllelic.go")
        intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
        TEST_DATA = intermediate.astype("object")

        _ = genomic_position_data._create_heatmap(
            TEST_DATA, min_values=1, height=600, width=600
        )

        mocked_go.Figure.assert_called_once()
        mocked_go.Heatmap.assert_called_once()

    def test_heatmap(self, set_up_genomic_position_data, mocker):
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.pyllelic.go")

        genomic_position_data.heatmap(min_values=1)

        mocked_go.Figure.assert_called_once()
        mocked_go.Heatmap.assert_called_once()

    def test_summarize_allelelic_data(self, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data
        EXPECTED = pd.DataFrame(
            {
                "cellLine": [],
                "position": [],
                "ad_stat": [],
                "p_crit": [],
                "diff": [],
                "raw": [],
            }
        )

        actual = genomic_position_data.summarize_allelic_data()
        # print(actual.to_markdown())

        pd.testing.assert_frame_equal(EXPECTED, actual)

    def test__create_methylation_diffs_bar_graph(
        self, mocker, set_up_genomic_position_data
    ):
        _, genomic_position_data = set_up_genomic_position_data
        mocked_go = mocker.patch("pyllelic.pyllelic.go")
        intermediate = pd.DataFrame.from_dict(
            {
                "cellLine": ["1", "2"],
                "position": ["1", "2"],
                "ad_stat": [0, 1],
                "p_crit": [1, 2],
                "diff": [1, 2],
                "raw": [1, 2],
            }
        )
        TEST_DATA = intermediate.astype("object")

        _ = genomic_position_data._create_methylation_diffs_bar_graph(TEST_DATA)

        mocked_go.Figure.assert_called_once()
        mocked_go.Bar.assert_called_once()

    # def test_sig_methylation_differences(self, set_up_genomic_position_data, mocker):
    #     _, genomic_position_data = set_up_genomic_position_data
    #     mocked_go = mocker.patch("pyllelic.pyllelic.go")

    #     genomic_position_data.sig_methylation_differences()

    #     mocked_go.Figure.assert_called_once()
    #     mocked_go.Bar.assert_called_once()

    def test_anderson_darling_test_with_values(self, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data
        TEST_LIST = pd.Series([1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)])
        actual = genomic_position_data._anderson_darling_test(TEST_LIST)
        EXPECTED = (
            False,
            0.9267475729317516,
            np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
        )
        assert EXPECTED[0] == actual[0]
        np.testing.assert_equal(EXPECTED[1], actual[1])
        assert (EXPECTED[2] == actual[2]).all()

    def test_anderson_darling_test_with_bad(self, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data
        TEST_LIST = [np.nan]
        actual = genomic_position_data._anderson_darling_test(TEST_LIST)
        EXPECTED = (False, np.nan, [np.nan])
        assert EXPECTED[0] == actual[0]
        np.testing.assert_equal(EXPECTED[1], actual[1])
        assert EXPECTED[2] == actual[2]

    def test_generate_ad_stats(self, set_up_genomic_position_data):
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

    def test_write_means_modes_diffs(self, mocker, set_up_genomic_position_data):
        _, genomic_position_data = set_up_genomic_position_data
        TEST_FILENAME = "output"
        mocker.patch.object(genomic_position_data.means, "to_excel")
        mocker.patch.object(genomic_position_data.modes, "to_excel")
        mocker.patch.object(genomic_position_data.diffs, "to_excel")

        genomic_position_data.write_means_modes_diffs(TEST_FILENAME)

        genomic_position_data.means.to_excel.assert_called()
        genomic_position_data.modes.to_excel.assert_called()
        genomic_position_data.diffs.to_excel.assert_called()

    # def test_save(self, mocker, set_up_genomic_position_data):
    #     _, genomic_position_data = set_up_genomic_position_data
    #     TEST_FILENAME = "output"
    #     mocked_excel_writer = mocker.patch("pyllelic.pyllelic.pd.ExcelWriter")

    #     genomic_position_data.save(TEST_FILENAME)

    #     mocked_excel_writer.ExcelWriter.assert_called()


# def test_return_individual_data():
#     """Check whether the expected and result DataFrames are identical."""
#     in_pos = ["1", "2", "3"]
#     in_cell = ["TEST1", "TEST2", "TEST3"]
#     result = pyllelic.return_individual_data(
#         dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
#     )
#     print(result.to_markdown())

#     intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA

#     expected = intermediate.astype("object")

#     pd.testing.assert_frame_equal(result, expected)


# def test_return_read_values():
#     """Check whether returned fractional methylation list is correct."""
#     in_pos = "1"
#     in_key = "TEST1"
#     in_min_reads = 1
#     in_min_sites = 1

#     expected = [
#         1.0,
#         1.0,
#         1.0,
#         2 / 3,
#         2 / 3,
#         2 / 3,
#     ]

#     result = pyllelic.return_read_values(
#         pos=in_pos,
#         key=in_key,
#         dict_of_dfs=SAMPLE_DICT_OF_DFS,
#         min_reads=in_min_reads,
#         min_sites=in_min_sites,
#     )

#     assert result == expected


# def test_get_str_values():
#     TEST_DATAFRAME = SAMPLE_DICT_OF_DFS.get("TEST1")
#     TEST_POSITION = "1"
#     EXPECTED = pd.Series(["111", "111", "111", "011", "011", "011"], name="1")
#     actual = pyllelic.get_str_values(TEST_DATAFRAME, TEST_POSITION)

#     pd.testing.assert_series_equal(EXPECTED, actual)


# def test_find_diffs():
#     """Check whether the expected and result DataFrames are identical."""
#     means = EXPECTED_INTERMEDIATE_MEANS
#     modes = EXPECTED_INTERMEDIATE_MODES

#     expected = EXPECTED_INTERMEDIATE_DIFFS

#     result = pyllelic.find_diffs(means, modes)

#     pd.testing.assert_frame_equal(result, expected)


# def test_truncate_diffs():
#     EXPECTED = EXPECTED_INTERMEDIATE_DIFFS.dropna(how="all")
#     actual = pyllelic.truncate_diffs(EXPECTED_INTERMEDIATE_DIFFS)
#     pd.testing.assert_frame_equal(EXPECTED, actual)
