#!/usr/bin/env python3
"""pytest unit tests for pyllelic."""

# Testing
import base64
import pytest
import unittest.mock as mock
import os
import tempfile
from contextlib import contextmanager

# import tempfile
from bam_input import (
    SAMPLE_BAM,
    SAMPLE_BAI,
    TEST_PROM_FILE,
    EXPECTED_BAM_OUTPUT_POSITIONS,
    EXPECTED_BAM_OUTPUT_VALUES,
    # SAMPLE_DICT_OF_DFS,
    # TEST_DF_SEQ_READS,
    # EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA,
    # EXPECTED_INTERMEDIATE_MODES,
    # EXPECTED_INTERMEDIATE_DIFFS,
    # EXPECTED_INTERMEDIATE_MEANS,
    EXPECTED_BAM_OUTPUT_GENOME_VALUES,
    EXPECTED_STACKED_BAM,
    EXPECTED_WRITE_DF_OUTPUT,
    INPUT_READ_FILE,
)

# Required libraries for test data
import pandas as pd

# import numpy as np
from pathlib import Path

# Module to test
import pyllelic.pyllelic as pyllelic


# Helper methods
@contextmanager
def tempinput(data):  # pragma: no cover
    """Helper for virtual files."""
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.write(data)
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)


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
    TEST_START = "1293000"
    TEST_END = "1296000"
    TEST_CHR = "5"
    TEST_OFFSET = 0
    pyllelic.set_up_env_variables(
        base_path=my_path,
        prom_file=prom_file,
        prom_start=TEST_START,
        prom_end=TEST_END,
        chrom=TEST_CHR,
        offset=TEST_OFFSET,
    )


# Tests
def test_set_up_env_variables():
    """Test setting environment variables with mock object."""
    TEST_BASE_PATH = Path().cwd()
    TEST_PROM_FILE = TEST_BASE_PATH / "test.txt"
    TEST_START = "1"
    TEST_END = "2"
    TEST_CHR = "5"
    TEST_OFFSET = 0
    EXPECTED_RESULTS = TEST_BASE_PATH / "results"
    pyllelic.set_up_env_variables(
        base_path=TEST_BASE_PATH,
        prom_file=TEST_PROM_FILE,
        prom_start=TEST_START,
        prom_end=TEST_END,
        chrom=TEST_CHR,
        offset=TEST_OFFSET,
    )

    assert TEST_BASE_PATH == pyllelic.config.base_directory
    assert TEST_PROM_FILE == pyllelic.config.promoter_file
    assert TEST_START == pyllelic.config.promoter_start
    assert TEST_END == pyllelic.config.promoter_end
    assert EXPECTED_RESULTS == pyllelic.config.results_directory


def test_make_list_of_bam_files():
    """Test making list of bam files, mocking config."""
    TEST_LIST = [
        Path("bad1.txt"),
        Path("bad2.bai"),
        Path("good1.bam"),
        Path("good2.bam"),
    ]
    EXPECTED = ["good1.bam", "good2.bam"]
    with mock.patch.object(pyllelic.Path, "iterdir") as mock_iterdir:
        mock_iterdir.return_value = TEST_LIST
        actual = pyllelic.make_list_of_bam_files()

    assert EXPECTED == actual


# Tests of main classes


class Test_BamOutput:
    """Class to test BamOutput object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_bam_output(self, tmp_path):
        p, fp_bam = setup_bam_files(tmp_path)
        setup_config(p)
        return pyllelic.BamOutput(
            sam_directory=fp_bam,
            genome_string=TEST_PROM_FILE,
            config=pyllelic.config,
        )

    def test_init(self, tmp_path, set_up_bam_output):
        bam_output = set_up_bam_output

        assert bam_output.name == str(tmp_path / "test" / "fh_test.bam")
        assert bam_output.values == EXPECTED_BAM_OUTPUT_VALUES
        pd.testing.assert_index_equal(
            bam_output.positions, EXPECTED_BAM_OUTPUT_POSITIONS
        )
        assert bam_output.genome_values == EXPECTED_BAM_OUTPUT_GENOME_VALUES

    def test_run_sam_and_extract_df(self, set_up_bam_output):
        bam_output = set_up_bam_output
        actual_positions = bam_output.run_sam_and_extract_df(Path(bam_output.name))
        pd.testing.assert_index_equal(actual_positions, EXPECTED_BAM_OUTPUT_POSITIONS)

    def test_write_bam_output(self, set_up_bam_output):
        bam_output = set_up_bam_output
        # print(f"Initial bam values:\n{bam_output.values}")
        bam_output.write_bam_output(bam_output.positions, EXPECTED_STACKED_BAM)
        # print(f"Revised bam values:\n{bam_output.values}")
        assert bam_output.values == EXPECTED_WRITE_DF_OUTPUT

    @mock.patch("pyllelic.pyllelic.pysam")
    def test_pysam_index(self, mock_pysam, set_up_bam_output):
        bam_output = set_up_bam_output
        TEST_PATH = Path().cwd()
        bam_output.pysam_index(TEST_PATH)
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
                    0: "1111",
                    1: "11111",
                    2: "11111",
                    3: "11111",
                    4: "11111",
                    5: "11111",
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
        EXPECTED = ["FAIL", "11"]
        actual = quma_output.process_raw_quma(TEST_QUMA_RESULT)
        assert EXPECTED == actual

    def test__pool_processing(self, set_up_quma_output):
        quma_output = set_up_quma_output

        EXPECTED_RESULTS = pd.DataFrame.from_dict(
            {
                "1293588": {
                    0: "1111",
                    1: "11111",
                    2: "11111",
                    3: "11111",
                    4: "11111",
                    5: "11111",
                }
            }
        )

        actual = quma_output._pool_processing()

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
            "genome\t0\tATCGTAGTCGA\t2\t2,8\n"
            + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
            + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
            + "2\tquery2\tATCGATAGCATT\tATCG-TAGT\tATCGATAGC\t"
            + "9\t1\t88.9\t1\t1\t0\t1\t100.0\t1\t1\t1\n"
        )
        # with mock.patch("builtins.open", new=my_open):
        actual = quma_output.access_quma(TEST_GSEQ, TEST_QSEQ)
        assert EXPECTED == actual


class Test_GenomicPositionData:
    """Class to test QumaOutput object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_genomic_position_data(self, tmp_path):
        p, _ = setup_bam_files(tmp_path)
        setup_config(p)
        INPUT_BAM_LIST = ["fh_test.bam"]
        return pyllelic.GenomicPositionData(
            config=pyllelic.config, files_set=INPUT_BAM_LIST
        )

    def test_init(self, set_up_genomic_position_data):
        genomic_position_data = set_up_genomic_position_data

        positions = []
        for each in EXPECTED_BAM_OUTPUT_POSITIONS:
            positions.append(each)

        EXPECTED_POSITIONS = sorted(set(positions))

        assert genomic_position_data.positions == EXPECTED_POSITIONS


# # Mock out main to ensure it just tests functionality of the main function
# @mock.patch("pyllelic.pyllelic.make_list_of_bam_files")
# @mock.patch("pyllelic.pyllelic.index_and_fetch")
# @mock.patch("pyllelic.pyllelic.genome_parsing")
# @mock.patch("pyllelic.pyllelic.extract_cell_types")
# @mock.patch("pyllelic.pyllelic.run_quma_and_compile_list_of_df")
# @mock.patch("pyllelic.pyllelic.process_means")
# @mock.patch("pyllelic.pyllelic.process_modes")
# @mock.patch("pyllelic.pyllelic.find_diffs")
# @mock.patch("pyllelic.pyllelic.write_means_modes_diffs")
# @mock.patch("pyllelic.pyllelic.sys.argv")
# def test_main(
#     mock_argv,
#     mock_writer,
#     mock_diffs,
#     mock_modes,
#     mock_means,
#     mock_run_quma,
#     mock_extract,
#     mock_genome_parse,
#     mock_index,
#     mock_bamlist,
# ):
#     """Test main module with a bunch of mocks."""
#     # Set up, patching all called functions
#     mock_bamlist.return_value = ["good1.bam", "good2.bam"]
#     mock_index.return_value = ["1", "2", "3"]
#     mock_genome_parse.return_value = None
#     mock_extract.return_value = ["TEST1", "TEST2", "TEST3"]
#     mock_run_quma.return_value = SAMPLE_DICT_OF_DFS
#     intermediate_means = EXPECTED_INTERMEDIATE_MEANS
#     mock_means.return_value = intermediate_means.astype("object")
#     intermediate_modes = EXPECTED_INTERMEDIATE_MODES
#     mock_modes.return_value = intermediate_modes.astype("object")
#     intermediate_diffs = EXPECTED_INTERMEDIATE_DIFFS
#     mock_diffs.return_value = intermediate_diffs.astype("object")
#     mock_writer.return_value = None
#     mock_argv.return_value = ["program", "output.xlsx"]

#     pyllelic.main()

#     # Asserts
#     mock_bamlist.assert_called_once()
#     mock_index.assert_called_once_with(["good1.bam", "good2.bam"])
#     mock_extract.assert_called_once_with(["good1.bam", "good2.bam"])
#     mock_run_quma.assert_called_once()
#     mock_means.assert_called_once()
#     mock_modes.assert_called_once()
#     mock_writer.assert_called_once()


# def test_extract_cell_types():
#     """Check correct response and exception to bad data."""
#     input_list = [
#         "fh_TEST1_TISSUE.TERT.bam",
#         "fh_TEST2_TISSUE.TERT.bam",
#         "fh_TEST3_TISSUE.TERT.bam",
#     ]
#     expected_result = ["TEST1", "TEST2", "TEST3"]
#     result = pyllelic.extract_cell_types(input_list)
#     assert expected_result == result

#     with pytest.raises(IndexError):
#         bad_input = ["text", "text2_", "_text3"]
#         pyllelic.extract_cell_types(bad_input)


# @mock.patch("pyllelic.pyllelic.read_df_of_quma_results")
# @mock.patch("pyllelic.pyllelic.quma_full_threaded")
# def test_run_quma_and_compile_list_of_df(mock_quma, mock_read):
#     TEST_CELL_TYPES = ["TEST1", "TEST2", "TEST3"]
#     TEST_FILENAME = "output.txt"
#     pyllelic.run_quma_and_compile_list_of_df(TEST_CELL_TYPES, TEST_FILENAME, True)
#     mock_quma.assert_called_once()
#     mock_read.assert_called_once()


# @mock.patch("pyllelic.pyllelic.pd")
# def test_read_df_of_quma_results(mock_pandas):
#     TEST_FILENAME = "output.xls"
#     pyllelic.read_df_of_quma_results(TEST_FILENAME)
#     mock_pandas.read_excel.assert_called_once()


# def test_process_means():
#     """Check whether the expected and result DataFrames are identical."""
#     in_pos = ["1", "2", "3"]
#     in_cell = ["TEST1", "TEST2", "TEST3"]
#     result = pyllelic.process_means(
#         dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
#     )

#     intermediate = EXPECTED_INTERMEDIATE_MEANS

#     expected = intermediate.astype("object")

#     pd.testing.assert_frame_equal(result, expected)


# def test_process_modes():
#     """Check whether the expected and result DataFrames are identical."""
#     in_pos = ["1", "2", "3"]
#     in_cell = ["TEST1", "TEST2", "TEST3"]
#     result = pyllelic.process_modes(
#         dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
#     )

#     intermediate = EXPECTED_INTERMEDIATE_MODES

#     expected = intermediate.astype("object")

#     pd.testing.assert_frame_equal(result, expected)


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


# # https://coderbook.com/@marcus/how-to-mock-and-unit-test-with-pandas/
# # https://stackoverflow.com/questions/65579240/unittest-mock-pandas-to-csv
# def test_write_means_modes_diffs():
#     TEST_DF = EXPECTED_INTERMEDIATE_MEANS
#     TEST_FILENAME = "output"
#     with mock.patch.object(TEST_DF, "to_excel") as to_excel_mock:
#         pyllelic.write_means_modes_diffs(TEST_DF, TEST_DF, TEST_DF, TEST_FILENAME)

#         to_excel_mock.assert_called()


# @mock.patch("pyllelic.pyllelic.go")
# def test_create_histogram(mock_go):
#     TEST_CELL_LINE = "TEST1"
#     TEST_POSITION = "1"
#     intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
#     TEST_DATA = intermediate.astype("object")

#     _ = pyllelic.create_histogram(TEST_DATA, TEST_CELL_LINE, TEST_POSITION)

#     mock_go.Figure.assert_called_once()
#     mock_go.Histogram.assert_called_once_with(
#         x=TEST_DATA.loc[TEST_CELL_LINE, TEST_POSITION],
#         xbins=dict(
#             start=-0.1,
#             end=1.1,
#             size=0.2,
#         ),
#     )


# @mock.patch("pyllelic.pyllelic.go")
# def test_histogram(mock_go):
#     TEST_CELL_LINE = "TEST1"
#     TEST_POSITION = "1"
#     intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
#     TEST_DATA = intermediate.astype("object")

#     pyllelic.histogram(TEST_DATA, TEST_CELL_LINE, TEST_POSITION)

#     mock_go.Figure.assert_called_once()
#     mock_go.Histogram.assert_called_once_with(
#         x=TEST_DATA.loc[TEST_CELL_LINE, TEST_POSITION],
#         xbins=dict(
#             start=-0.1,
#             end=1.1,
#             size=0.2,
#         ),
#     )


# def test_anderson_darling_test_with_values():
#     TEST_LIST = pd.Series([1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)])
#     actual = pyllelic.anderson_darling_test(TEST_LIST)
#     EXPECTED = (
#         False,
#         0.9267475729317516,
#         np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
#     )
#     assert EXPECTED[0] == actual[0]
#     np.testing.assert_equal(EXPECTED[1], actual[1])
#     assert (EXPECTED[2] == actual[2]).all()


# def test_anderson_darling_test_with_bad():
#     TEST_LIST = [np.nan]
#     actual = pyllelic.anderson_darling_test(TEST_LIST)
#     EXPECTED = (False, np.nan, [np.nan])
#     assert EXPECTED[0] == actual[0]
#     np.testing.assert_equal(EXPECTED[1], actual[1])
#     assert EXPECTED[2] == actual[2]


# def test_generate_ad_stats():
#     TEST_DF = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
#     EXPECTED = pd.DataFrame.from_dict(
#         {
#             "TEST1": [
#                 (
#                     False,
#                     0.9267475729317516,
#                     np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
#                 ),
#                 (
#                     False,
#                     0.9267475729317516,
#                     np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
#                 ),
#                 (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
#             ],
#             "TEST2": [
#                 (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
#                 (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
#                 (False, np.nan, np.array([0.592, 0.675, 0.809, 0.944, 1.123])),
#             ],
#             "TEST3": [
#                 (False, np.nan, [np.nan]),
#                 (False, np.nan, [np.nan]),
#                 (
#                     False,
#                     0.7995458667216608,
#                     np.array([0.72, 0.82, 0.984, 1.148, 1.365]),
#                 ),
#             ],
#         },
#         orient="index",
#         columns=["1", "2", "3"],
#     )
#     actual = pyllelic.generate_ad_stats(TEST_DF)
#     EXPECTED_TEST2_1 = EXPECTED.loc["TEST2", "1"]
#     np.testing.assert_equal(EXPECTED_TEST2_1, actual.loc["TEST2", "1"])


# def test_summarize_allelelic_data():
#     TEST_INDIV = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
#     TEST_DIFFS = EXPECTED_INTERMEDIATE_DIFFS
#     EXPECTED = pd.DataFrame(
#         {
#             "cellLine": [],
#             "position": [],
#             "ad_stat": [],
#             "p_crit": [],
#             "diff": [],
#             "raw": [],
#         }
#     )

#     actual = pyllelic.summarize_allelic_data(TEST_INDIV, TEST_DIFFS)
#     print(actual.to_markdown())

#     pd.testing.assert_frame_equal(EXPECTED, actual)
