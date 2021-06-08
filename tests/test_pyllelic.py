#!/usr/bin/env python3
"""pytest unit tests for pyllelic."""

# Testing
import pytest
import unittest.mock as mock
import os
import tempfile
from contextlib import contextmanager

# import tempfile

# Required libraries for test data
import pandas as pd
import numpy as np
from pathlib import Path

# Module to test
import pyllelic.pyllelic as pyllelic


# Test data
SAMPLE_DICT_OF_DFS = {
    "TEST1": pd.DataFrame(
        {
            "1": ["111", "111", "111", "011", "011", "011"],
            "2": ["111", "111", "111", "011", "011", "011"],
            "3": ["111", "111", "111", "111", "111", "111"],
        }
    ),
    "TEST2": pd.DataFrame(
        {
            "1": ["111", "111", "111", "111", "111", "111"],
            "2": ["111", "111", "111", "111", "111", "111"],
            "3": ["111", "111", "111", "111", "111", "111"],
        }
    ),
    "TEST3": pd.DataFrame(
        {
            "1": ["FAIL", "FAIL", "FAIL", "FAIL", "FAIL", "FAIL"],
            "2": ["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"],
            "3": ["111", "111", "111", "011", "011", "NaN"],
        }
    ),
}

TEST_DF_SEQ_READS = pd.DataFrame.from_dict(
    {
        "1": "ATGCGTACGTA",
        "2": "ATGCTAGGCCC",
    },
    orient="index",
    columns=["sequence"],
).stack()

EXPECTED_INTERMEDIATE_MEANS = pd.DataFrame.from_dict(
    {
        "TEST1": [np.float64(5 / 6), np.float64(5 / 6), np.float64(1.0)],
        "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
        "TEST3": [np.nan, np.nan, np.float64(13 / 15)],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_INTERMEDIATE_MODES = pd.DataFrame.from_dict(
    {
        "TEST1": [np.float64(2 / 3), np.float64(2 / 3), np.float64(1.0)],
        "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
        "TEST3": [np.nan, np.nan, np.float64(1.0)],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_INTERMEDIATE_DIFFS = pd.DataFrame.from_dict(
    {
        "TEST1": [np.float64(1 / 6), np.float64(1 / 6), np.float64(0.0)],
        "TEST2": [np.float64(0.0), np.float64(0.0), np.float64(0.0)],
        "TEST3": [np.nan, np.nan, np.float64(-2 / 15)],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA = pd.DataFrame.from_dict(
    {
        "TEST1": [
            [1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)],
            [1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ],
        "TEST2": [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ],
        "TEST3": [np.nan, np.nan, [1.0, 1.0, 1.0, (2 / 3), (2 / 3)]],
    },
    orient="index",
    columns=["1", "2", "3"],
)


# Helper method
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


# Mock out main to ensure it just tests functionality of the main function
@mock.patch("pyllelic.pyllelic.make_list_of_bam_files")
@mock.patch("pyllelic.pyllelic.index_and_fetch")
@mock.patch("pyllelic.pyllelic.genome_parsing")
@mock.patch("pyllelic.pyllelic.extract_cell_types")
@mock.patch("pyllelic.pyllelic.run_quma_and_compile_list_of_df")
@mock.patch("pyllelic.pyllelic.process_means")
@mock.patch("pyllelic.pyllelic.process_modes")
@mock.patch("pyllelic.pyllelic.find_diffs")
@mock.patch("pyllelic.pyllelic.write_means_modes_diffs")
@mock.patch("pyllelic.pyllelic.sys.argv")
def test_main(
    mock_argv,
    mock_writer,
    mock_diffs,
    mock_modes,
    mock_means,
    mock_run_quma,
    mock_extract,
    mock_genome_parse,
    mock_index,
    mock_bamlist,
):
    """Test main module with a bunch of mocks."""
    # Set up, patching all called functions
    mock_bamlist.return_value = ["good1.bam", "good2.bam"]
    mock_index.return_value = ["1", "2", "3"]
    mock_genome_parse.return_value = None
    mock_extract.return_value = ["TEST1", "TEST2", "TEST3"]
    mock_run_quma.return_value = SAMPLE_DICT_OF_DFS
    intermediate_means = EXPECTED_INTERMEDIATE_MEANS
    mock_means.return_value = intermediate_means.astype("object")
    intermediate_modes = EXPECTED_INTERMEDIATE_MODES
    mock_modes.return_value = intermediate_modes.astype("object")
    intermediate_diffs = EXPECTED_INTERMEDIATE_DIFFS
    mock_diffs.return_value = intermediate_diffs.astype("object")
    mock_writer.return_value = None
    mock_argv.return_value = ["program", "output.xlsx"]

    pyllelic.main()

    # Asserts
    mock_bamlist.assert_called_once()
    mock_index.assert_called_once_with(["good1.bam", "good2.bam"])
    mock_extract.assert_called_once_with(["good1.bam", "good2.bam"])
    mock_run_quma.assert_called_once()
    mock_means.assert_called_once()
    mock_modes.assert_called_once()
    mock_writer.assert_called_once()


def test_genome_range():
    """Check if correct genome string is returned."""
    gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    result = pyllelic.genome_range(position=2, genome_string=gen_str, offset=40)
    assert result == gen_str[8:37]
    assert isinstance(result, str)


def test_genome_range_no_offset():
    """Check if correct genome string is returned with no offset given."""
    gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    pyllelic.config.offset = 40
    result = pyllelic.genome_range(position=2, genome_string=gen_str)
    assert result == gen_str[8:37]
    assert isinstance(result, str)


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


@mock.patch("pyllelic.pyllelic.run_sam_and_extract_df")
def test_index_and_fetch(mock_run_sam):
    TEST_FILES = ["TEST1", "TEST2", "TEST3"]
    TEST_PROCESS = False
    with tempfile.TemporaryDirectory() as tmpdirname:
        pyllelic.config.base_directory = Path(tmpdirname)
        mock_run_sam.return_value = pd.Index(["1", "2", "3"])
        pyllelic.index_and_fetch(TEST_FILES, TEST_PROCESS)

        mock_run_sam.assert_called()


# @mock.patch("pyllelic.pyllelic.pysam")
# def test_run_sam_and_extract_df_no_process(mock_pysam):
#     """Check if a samfile will be properly aligned and read."""
#     TEST_SAM_FILE = Path("TEST1.bam")
#     TEST_SAM_ALIGNMENT =
#     EXPECTED = pd.Index([])
#     with mock.patch.object(pyllelic.Path, "exists") as mock_exists:
#         mock_exists.return_value = True
#         mock_pysam.AlignmentFile.return_value = TEST_SAM_ALIGNMENT
#         result = pyllelic.run_sam_and_extract_df(TEST_SAM_FILE, process=False)

#     mock_pysam.AlignmentFile.assert_called_once()

#     pd.testing.assert_index_equal(result, expected)


def test_run_sam_and_extract_df_with_process():
    """Check if a samfile will be properly aligned and read."""


@mock.patch("pyllelic.pyllelic.write_individual_bam_file")
def test_write_bam_output_files(mock_writer):
    """Test writing bam files with a mock open."""
    TEST_SAMS = Path("fh_TEST1/")
    TEST_POSITIONS = ["1", "2"]
    TEST_DF = TEST_DF_SEQ_READS
    with tempfile.TemporaryDirectory() as tmpdirname:
        pyllelic.config.base_directory = Path(tmpdirname)
        pyllelic.write_bam_output_files(TEST_SAMS, TEST_POSITIONS, TEST_DF)

    mock_writer.assert_called()


def test_write_individual_bam_file(mocker):
    """Check if bam outputs would be correctly written."""
    open_mock = mocker.mock_open()
    test_position = "1"
    test_contents = [">read0", "ATGCATGCATGCATGC"]
    test_sam = "fh_TEST1_CELL.TERT.bam"
    expected_directory = pyllelic.config.base_directory.joinpath(
        "bam_output", test_sam, "1.txt"
    )

    with mock.patch("builtins.open", open_mock, create=True) as m:
        pyllelic.write_individual_bam_file(
            sam_name=test_sam, filename=test_position, file_contents=test_contents
        )

    open_mock.assert_called_with(expected_directory, "w")
    handle = m()
    handle.write.assert_called_with("ATGCATGCATGCATGC\n")


@mock.patch("pyllelic.pyllelic.pysam")
def test_pysam_index(mock_pysam):
    TEST_PATH = Path().cwd()
    pyllelic.pysam_index(TEST_PATH)
    mock_pysam.index.assert_called_once_with(os.fspath(TEST_PATH))


@mock.patch("pyllelic.pyllelic._process_genome_parsing")
@mock.patch.object(pyllelic.Path, "is_dir")
@mock.patch.object(pyllelic.Path, "iterdir")
def test_genome_parsing(mock_iterdir, mock_is_dir, mock_process):
    with tempfile.NamedTemporaryFile() as my_file:
        my_file.write(b"ATGCTCTAGCTCGCTAGATCGCTCGATCGTAGCTAGCTA")
        my_file.seek(0)
        pyllelic.config.promoter_file = my_file.name
        mock_iterdir.return_value = [Path("TEST1/"), Path("TEST2/")]
        mock_is_dir.return_value = True
        pyllelic.genome_parsing()

        mock_process.assert_called()


@mock.patch("pyllelic.pyllelic.genome_range")
@mock.patch.object(pyllelic.os, "listdir")
def test__process_genome_parsing(mock_listdir, mock_range):
    with tempfile.TemporaryDirectory() as tempdir:
        TEST_GENOME_STRING = "ATGCTCTAGCTCGCTAGATCGCTCGATCGTAGCTAGCTA"
        mock_listdir.return_value = ["1.txt", "2.txt", "g_1.txt"]
        mock_range.return_value = "ACGACATGA"
        pyllelic._process_genome_parsing(Path(tempdir), TEST_GENOME_STRING)

        mock_range.assert_called()


def test_access_quma():
    TEST_DIRECTORY = "Test"
    TEST_GSEQ_NAME = "genome.txt"
    TEST_GSEQ = ">genome\nATCGTAGTCGA"
    TEST_QSEQ_NAME = "query.txt"
    TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"

    # From https://stackoverflow.com/questions/26783678/
    #      python-mock-builtin-open-in-a-class-using-two-different-files
    def my_open(filename, _):
        if filename == f"{TEST_DIRECTORY}/{TEST_GSEQ_NAME}":
            content = TEST_GSEQ
        elif filename == f"{TEST_DIRECTORY}/{TEST_QSEQ_NAME}":
            content = TEST_QSEQ
        else:  # pragma: no cover
            raise FileNotFoundError(filename)
        file_object = mock.mock_open(read_data=content).return_value
        file_object.__iter__.return_value = content.splitlines(True)
        return file_object

    EXPECTED = (
        "genome\t0\tATCGTAGTCGA\t2\t2,8\n"
        + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
        + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
        + "2\tquery2\tATCGATAGCATT\tATCG-TAGT\tATCGATAGC\t"
        + "9\t1\t88.9\t1\t1\t0\t1\t100.0\t1\t1\t1\n"
    )
    with mock.patch("builtins.open", new=my_open):
        actual = pyllelic.access_quma(TEST_DIRECTORY, TEST_GSEQ_NAME, TEST_QSEQ_NAME)
    assert EXPECTED == actual


@mock.patch("pyllelic.pyllelic.signal.signal")
def test__init_worker(mock_signal):
    pyllelic._init_worker()
    mock_signal.assert_called_once()


def test__thread_worker():
    TEST_FOLDER = "Test"
    TEST_READ_NAME = "1295094"
    TEST_GSEQ_NAME = f"g_{TEST_READ_NAME}.txt"
    TEST_GSEQ = ">genome\nATCGTAGTCGA"
    TEST_QSEQ_NAME = f"{TEST_READ_NAME}.txt"
    TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"

    EXPECTED: pd.DataFrame = pd.DataFrame({"1295094": ["1", "11"]})

    # From https://stackoverflow.com/questions/26783678/
    #      python-mock-builtin-open-in-a-class-using-two-different-files
    def my_open(filename, _):
        if filename == f"{TEST_FOLDER}/{TEST_GSEQ_NAME}":
            content = TEST_GSEQ
        elif filename == f"{TEST_FOLDER}/{TEST_QSEQ_NAME}":
            content = TEST_QSEQ
        else:  # pragma: no cover
            raise FileNotFoundError(filename)
        file_object = mock.mock_open(read_data=content).return_value
        file_object.__iter__.return_value = content.splitlines(True)
        return file_object

    with mock.patch("builtins.open", new=my_open):
        actual = pyllelic._thread_worker(TEST_FOLDER, TEST_READ_NAME)

    pd.testing.assert_frame_equal(actual, EXPECTED)


@mock.patch("pyllelic.pyllelic._pool_processing")
@mock.patch.object(pyllelic.os, "listdir")
@mock.patch.object(pyllelic.Path, "iterdir")
def test_quma_full_threaded(mock_iterdir, mock_listdir, mock_pool):
    TEST_TYPES = ["TEST1", "TEST2", "TEST3"]
    TEST_FILENAME = "output.xls"
    TEST_ITERDIR = [
        Path("bad1.txt"),
        Path("bad2.bai"),
        Path("fh_TEST1/"),
        Path("fh_TEST2/"),
    ]
    with tempfile.TemporaryDirectory() as tmpdirname:
        pyllelic.config.base_directory = Path(tmpdirname)
        mock_iterdir.return_value = TEST_ITERDIR
        mock_listdir.return_value = ["g_1.txt", "g_2.txt", "1.txt", "2.txt"]
        mock_pool.return_value = pd.DataFrame({"1295094": ["1A", "11"]})
        pyllelic.quma_full_threaded(TEST_TYPES, TEST_FILENAME)

        assert Path(tmpdirname).joinpath("output.xls").exists()


class MockPoolApplyResult:
    def __init__(self, func, args):
        self._func = func
        self._args = args

    def get(self, timeout=0):
        return self._func(*self._args)


def _mock_apply_async(
    self, func, args=(), kwds=None, callback=None, error_callback=None
):
    return MockPoolApplyResult(func, args)


@pytest.fixture(autouse=True)
def mock_pool_apply_async(monkeypatch):
    monkeypatch.setattr("multiprocessing.pool.Pool.apply_async", _mock_apply_async)


@mock.patch("pyllelic.pyllelic._thread_worker")
def test__pool_processing(mock_thread, mock_pool_apply_async):
    TEST_READ_FILES = ["1.txt", "2.txt"]
    EXPECTED_THREAD = pd.DataFrame({"1295094": ["1", "11"]})
    EXPECTED_RESULTS = pd.concat([EXPECTED_THREAD, EXPECTED_THREAD], axis=1)
    with tempfile.TemporaryDirectory() as tmpdirname:
        TEST_FOLDER = Path(tmpdirname)
        mock_thread.return_value = EXPECTED_THREAD
        actual = pyllelic._pool_processing(TEST_READ_FILES, TEST_FOLDER)

    pd.testing.assert_frame_equal(actual, EXPECTED_RESULTS)


def test_process_raw_quma():
    TEST_QUMA_RESULT = (
        "genome\t0\tATCGTAGTCGA\t2\t2,8\n"
        + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
        + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
        + "2\tquery2\tATCGATAGCATT\tATCGATAGCATT\tATCG-TAGTCGA\t"
        + "12\t5\t58.3\t1\t1\t0\t1\t100.0\t1A\t1\t1\n"
    )
    EXPECTED = ["FAIL", "11"]
    actual = pyllelic.process_raw_quma(TEST_QUMA_RESULT)
    assert EXPECTED == actual


def test_extract_cell_types():
    """Check correct response and exception to bad data."""
    input_list = [
        "fh_TEST1_TISSUE.TERT.bam",
        "fh_TEST2_TISSUE.TERT.bam",
        "fh_TEST3_TISSUE.TERT.bam",
    ]
    expected_result = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.extract_cell_types(input_list)
    assert expected_result == result

    with pytest.raises(IndexError):
        bad_input = ["text", "text2_", "_text3"]
        pyllelic.extract_cell_types(bad_input)


@mock.patch("pyllelic.pyllelic.read_df_of_quma_results")
@mock.patch("pyllelic.pyllelic.quma_full_threaded")
def test_run_quma_and_compile_list_of_df(mock_quma, mock_read):
    TEST_CELL_TYPES = ["TEST1", "TEST2", "TEST3"]
    TEST_FILENAME = "output.txt"
    pyllelic.run_quma_and_compile_list_of_df(TEST_CELL_TYPES, TEST_FILENAME, True)
    mock_quma.assert_called_once()
    mock_read.assert_called_once()


@mock.patch("pyllelic.pyllelic.pd")
def test_read_df_of_quma_results(mock_pandas):
    TEST_FILENAME = "output.xls"
    pyllelic.read_df_of_quma_results(TEST_FILENAME)
    mock_pandas.read_excel.assert_called_once()


def test_process_means():
    """Check whether the expected and result DataFrames are identical."""
    in_pos = ["1", "2", "3"]
    in_cell = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.process_means(
        dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    )

    intermediate = EXPECTED_INTERMEDIATE_MEANS

    expected = intermediate.astype("object")

    pd.testing.assert_frame_equal(result, expected)


def test_process_modes():
    """Check whether the expected and result DataFrames are identical."""
    in_pos = ["1", "2", "3"]
    in_cell = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.process_modes(
        dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    )

    intermediate = EXPECTED_INTERMEDIATE_MODES

    expected = intermediate.astype("object")

    pd.testing.assert_frame_equal(result, expected)


def test_return_individual_data():
    """Check whether the expected and result DataFrames are identical."""
    in_pos = ["1", "2", "3"]
    in_cell = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.return_individual_data(
        dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    )
    print(result.to_markdown())

    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA

    expected = intermediate.astype("object")

    pd.testing.assert_frame_equal(result, expected)


def test_return_read_values():
    """Check whether returned fractional methylation list is correct."""
    in_pos = "1"
    in_key = "TEST1"
    in_min_reads = 1
    in_min_sites = 1

    expected = [
        1.0,
        1.0,
        1.0,
        2 / 3,
        2 / 3,
        2 / 3,
    ]

    result = pyllelic.return_read_values(
        pos=in_pos,
        key=in_key,
        dict_of_dfs=SAMPLE_DICT_OF_DFS,
        min_reads=in_min_reads,
        min_sites=in_min_sites,
    )

    assert result == expected


def test_get_str_values():
    TEST_DATAFRAME = SAMPLE_DICT_OF_DFS.get("TEST1")
    TEST_POSITION = "1"
    EXPECTED = pd.Series(["111", "111", "111", "011", "011", "011"], name="1")
    actual = pyllelic.get_str_values(TEST_DATAFRAME, TEST_POSITION)

    pd.testing.assert_series_equal(EXPECTED, actual)


def test_find_diffs():
    """Check whether the expected and result DataFrames are identical."""
    means = EXPECTED_INTERMEDIATE_MEANS
    modes = EXPECTED_INTERMEDIATE_MODES

    expected = EXPECTED_INTERMEDIATE_DIFFS

    result = pyllelic.find_diffs(means, modes)

    pd.testing.assert_frame_equal(result, expected)


def test_truncate_diffs():
    EXPECTED = EXPECTED_INTERMEDIATE_DIFFS.dropna(how="all")
    actual = pyllelic.truncate_diffs(EXPECTED_INTERMEDIATE_DIFFS)
    pd.testing.assert_frame_equal(EXPECTED, actual)


# https://coderbook.com/@marcus/how-to-mock-and-unit-test-with-pandas/
# https://stackoverflow.com/questions/65579240/unittest-mock-pandas-to-csv
def test_write_means_modes_diffs():
    TEST_DF = EXPECTED_INTERMEDIATE_MEANS
    TEST_FILENAME = "output"
    with mock.patch.object(TEST_DF, "to_excel") as to_excel_mock:
        pyllelic.write_means_modes_diffs(TEST_DF, TEST_DF, TEST_DF, TEST_FILENAME)

        to_excel_mock.assert_called()


@mock.patch("pyllelic.pyllelic.go")
def test_create_histogram(mock_go):
    TEST_CELL_LINE = "TEST1"
    TEST_POSITION = "1"
    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    _ = pyllelic.create_histogram(TEST_DATA, TEST_CELL_LINE, TEST_POSITION)

    mock_go.Figure.assert_called_once()
    mock_go.Histogram.assert_called_once_with(
        x=TEST_DATA.loc[TEST_CELL_LINE, TEST_POSITION],
        xbins=dict(
            start=-0.1,
            end=1.1,
            size=0.2,
        ),
    )


@mock.patch("pyllelic.pyllelic.go")
def test_histogram(mock_go):
    TEST_CELL_LINE = "TEST1"
    TEST_POSITION = "1"
    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    pyllelic.histogram(TEST_DATA, TEST_CELL_LINE, TEST_POSITION)

    mock_go.Figure.assert_called_once()
    mock_go.Histogram.assert_called_once_with(
        x=TEST_DATA.loc[TEST_CELL_LINE, TEST_POSITION],
        xbins=dict(
            start=-0.1,
            end=1.1,
            size=0.2,
        ),
    )


def test_anderson_darling_test_with_values():
    TEST_LIST = pd.Series([1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)])
    actual = pyllelic.anderson_darling_test(TEST_LIST)
    EXPECTED = (
        False,
        0.9267475729317516,
        np.array([0.592, 0.675, 0.809, 0.944, 1.123]),
    )
    assert EXPECTED[0] == actual[0]
    np.testing.assert_equal(EXPECTED[1], actual[1])
    assert (EXPECTED[2] == actual[2]).all()


def test_anderson_darling_test_with_bad():
    TEST_LIST = [np.nan]
    actual = pyllelic.anderson_darling_test(TEST_LIST)
    EXPECTED = (False, np.nan, [np.nan])
    assert EXPECTED[0] == actual[0]
    np.testing.assert_equal(EXPECTED[1], actual[1])
    assert EXPECTED[2] == actual[2]


def test_generate_ad_stats():
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
    actual = pyllelic.generate_ad_stats(TEST_DF)
    EXPECTED_TEST2_1 = EXPECTED.loc["TEST2", "1"]
    np.testing.assert_equal(EXPECTED_TEST2_1, actual.loc["TEST2", "1"])


def test_summarize_allelelic_data():
    TEST_INDIV = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DIFFS = EXPECTED_INTERMEDIATE_DIFFS
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

    actual = pyllelic.summarize_allelic_data(TEST_INDIV, TEST_DIFFS)
    print(actual.to_markdown())

    pd.testing.assert_frame_equal(EXPECTED, actual)
