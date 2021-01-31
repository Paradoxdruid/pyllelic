#!/usr/bin/env python3
"""pytest unit tests for pyllelic."""

# Testing
import pytest
import unittest.mock as mock

# import tempfile

# Required libraries for test data
import pandas as pd
import numpy as np
from pathlib import Path

# Module to test
import pyllelic


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


# Tests
def test_set_up_env_variables(mocker):
    """Test setting environment variables with mock object."""
    # """Check if config variables would be properly set."""
    # sys.modules["pyllelic.config"] = mocker.MagicMock()
    # pyllelic.config.base_directory = mock.MagicMock()
    mocker.patch.object(pyllelic, "config", autospec=True)
    # mocker.patch.object(pyllelic.config, "base_directory")
    # mocker.patch.object(pyllelic.config, "promoter_file")
    # mocker.patch.object(pyllelic.config, "results_directory")
    # mocker.patch.object(pyllelic.config, "bam_directory")
    # mocker.patch.object(pyllelic.config, "analysis_directory")
    # mocker.patch.object(pyllelic.config, "promoter_start")
    # mocker.patch.object(pyllelic.config, "promoter_end")
    # mocker.patch.object(pyllelic.config, "chromosome")

    pyllelic.set_up_env_variables(
        base_path=".", prom_file="test.txt", prom_start="1", prom_end="2", chrom="5"
    )

    # pyllelic.config.base_directory.assert_called()  # _with(Path("."))
    # assert pyllelic.config.base_directory == Path(".")
    assert pyllelic.config.base_directory == Path(".")


def test_main(mocker):
    """Test main module with a bunch of mocks."""
    # Set up, patching all called functions
    # mocker.patch("pyllelic.sys.argv", return_value=["program", "output.xlsx"])

    # expected_bam_files_return = [
    #     "fh_TEST1_TISSUE.TERT.bam",
    #     "fh_TEST2_TISSUE.TERT.bam",
    #     "fh_TEST3_TISSUE.TERT.BAM",
    # ]
    # mocker.patch(
    #     "pyllelic.make_list_of_bam_files",
    #     return_value=expected_bam_files_return,
    # )
    # expected_positions = ["1", "2", "3"]
    # mocker.patch("pyllelic.index_and_fetch", return_value=expected_positions)
    # mocker.patch("pyllelic.genome_parsing", return_value=None)
    # expected_cell_types = ["TEST1", "TEST2"]
    # mocker.patch("pyllelic.extract_cell_types", return_value=expected_cell_types)
    # mocker.patch(
    #     "pyllelic.run_quma_and_compile_list_of_df", return_value=SAMPLE_DICT_OF_DFS
    # )

    # means_intermediate = pd.DataFrame.from_dict(
    #     {
    #         "TEST1": [np.float64(5 / 6), np.float64(5 / 6), np.float64(1.0)],
    #         "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
    #         "TEST3": [np.nan, np.nan, np.float64(13 / 15)],
    #     },
    #     orient="index",
    #     columns=["1", "2", "3"],
    # )

    # means_expected = means_intermediate.astype("object")
    # mocker.patch("pyllelic.process_means", return_value=means_expected)

    # modes_intermediate = pd.DataFrame.from_dict(
    #     {
    #         "TEST1": [np.float64(2 / 3), np.float64(2 / 3), np.float64(1.0)],
    #         "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
    #         "TEST3": [np.nan, np.nan, np.float64(1.0)],
    #     },
    #     orient="index",
    #     columns=["1", "2", "3"],
    # )

    # modes_expected = modes_intermediate.astype("object")
    # mocker.patch("pyllelic.process_modes", return_value=modes_expected)

    # diffs_expected = pd.DataFrame.from_dict(
    #     {
    #         "TEST1": [np.float64(1 / 6), np.float64(1 / 6), np.float64(0.0)],
    #         "TEST2": [np.float64(0.0), np.float64(0.0), np.float64(0.0)],
    #         "TEST3": [np.nan, np.nan, np.float64(-2 / 15)],
    #     },
    #     orient="index",
    #     columns=["1", "2", "3"],
    # )
    # mocker.patch("pyllelic.find_diffs", return_value=diffs_expected)
    # mocker.patch("pyllelic.write_means_modes_diffs")

    # # Run it
    # pyllelic.main()

    # # Assertions at each step
    # pyllelic.make_list_of_bam_files.assert_called_once_with()
    # pyllelic.index_and_fetch.assert_called_once_with(expected_bam_files_return)
    # pyllelic.genome_parsing.assert_called_once_with()
    # pyllelic.extract_cell_types.assert_called_once_with(expected_bam_files_return)
    # pyllelic.run_quma_and_compile_list_of_df.assert_called_once_with(
    #     expected_cell_types, "output.xlsx"
    # )
    # pyllelic.process_means.assert_called_once_with(
    #     SAMPLE_DICT_OF_DFS, expected_positions, expected_cell_types
    # )
    # pyllelic.process_modes.assert_called_once_with(
    #     SAMPLE_DICT_OF_DFS, expected_positions, expected_cell_types
    # )
    # pyllelic.find_diffs.assert_called_once_with(means_expected, modes_expected)
    # pyllelic.write_means_modes_diffs.assert_called_once_with(
    #     means_expected, modes_expected, diffs_expected, "output.xlsx"
    # )


def test_genome_range():
    """Check if correct genome string is returned."""
    gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    result = pyllelic.genome_range(position=2, genome_string=gen_str, offset=40)
    assert result == gen_str[8:37]
    assert isinstance(result, str)


def test_make_list_of_bam_files():
    """Test making list of bam files, mocking config."""
    # with mock.patch("pyllelic.config.analysis_directory") as mock_dir:
    #     mock_dir.__get__ = mock.Mock(return_value=Path())


def test_index_and_fetch():
    pass


def test_run_sam_and_extract_df(mocker):
    """Check if a samfile will be properly aligned and read."""
    # TEST_SAM_FILE = Path("TEST1")
    # mocker.patch.object(pyllelic.pysam, "AlignmentFile", autospec=True)
    # mocker.patch.object(pyllelic.config, "base_directory", autospec=True)
    # expected = pd.Index([], name="positions")

    # result = pyllelic.run_sam_and_extract_df(TEST_SAM_FILE)

    # pd.testing.assert_index_equal(result, expected)


def test_write_bam_output_files(mocker):
    """Test writing bam files with a mock open."""
    # with tempfile.TemporaryDirectory() as tmpdirname:
    #     mocker.patch.object(pyllelic.config, "base_directory", Path(tmpdirname))


def test_write_individual_bamfile(mocker):
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


def test_samtools_index():
    pass


def test_genome_parsing():
    pass


def test_run_quma():
    pass


def test_quma_full():
    pass


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


def test_run_quma_and_compile_list_of_df():
    pass


def test_read_df_of_quma_results():
    pass


def test_process_means():
    """Check whether the expected and result DataFrames are identical."""
    in_pos = ["1", "2", "3"]
    in_cell = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.process_means(
        dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    )

    intermediate = pd.DataFrame.from_dict(
        {
            "TEST1": [np.float64(5 / 6), np.float64(5 / 6), np.float64(1.0)],
            "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
            "TEST3": [np.nan, np.nan, np.float64(13 / 15)],
        },
        orient="index",
        columns=["1", "2", "3"],
    )

    expected = intermediate.astype("object")

    pd.testing.assert_frame_equal(result, expected)


def test_process_modes():
    """Check whether the expected and result DataFrames are identical."""
    in_pos = ["1", "2", "3"]
    in_cell = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.process_modes(
        dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    )

    intermediate = pd.DataFrame.from_dict(
        {
            "TEST1": [np.float64(2 / 3), np.float64(2 / 3), np.float64(1.0)],
            "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
            "TEST3": [np.nan, np.nan, np.float64(1.0)],
        },
        orient="index",
        columns=["1", "2", "3"],
    )

    expected = intermediate.astype("object")

    pd.testing.assert_frame_equal(result, expected)


def test_return_individual_data():
    """Check whether the expected and result DataFrames are identical."""
    #     in_pos = ["1", "2", "3"]
    #     in_cell = ["TEST1", "TEST2", "TEST3"]
    #     result = pyllelic.return_individual_data(
    #         dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    #     )

    #     intermediate = pd.DataFrame.from_dict(
    #         {
    #             "TEST1": [
    #                 [1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)],
    #                 [1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)],
    #                 [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #             ],
    #             "TEST2": [
    #                 [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #                 [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #                 [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #             ],
    #             "TEST3": [
    #                 [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
    #                 [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
    #                 [1.0, 1.0, 1.0, (2 / 3), (2 / 3), np.nan],
    #             ],
    #         },
    #         orient="index",
    #         columns=["1", "2", "3"],
    #     )

    #     expected = intermediate.astype("object")

    #     pd.testing.assert_frame_equal(result, expected)


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


def test_find_diffs():
    """Check whether the expected and result DataFrames are identical."""
    means = pd.DataFrame.from_dict(
        {
            "TEST1": [np.float64(5 / 6), np.float64(5 / 6), np.float64(1.0)],
            "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
            "TEST3": [np.nan, np.nan, np.float64(13 / 15)],
        },
        orient="index",
        columns=["1", "2", "3"],
    )

    modes = pd.DataFrame.from_dict(
        {
            "TEST1": [np.float64(2 / 3), np.float64(2 / 3), np.float64(1.0)],
            "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
            "TEST3": [np.nan, np.nan, np.float64(1.0)],
        },
        orient="index",
        columns=["1", "2", "3"],
    )

    expected = pd.DataFrame.from_dict(
        {
            "TEST1": [np.float64(1 / 6), np.float64(1 / 6), np.float64(0.0)],
            "TEST2": [np.float64(0.0), np.float64(0.0), np.float64(0.0)],
            "TEST3": [np.nan, np.nan, np.float64(-2 / 15)],
        },
        orient="index",
        columns=["1", "2", "3"],
    )

    result = pyllelic.find_diffs(means, modes)

    pd.testing.assert_frame_equal(result, expected)


def test_write_means_modes_diffs():
    pass


def test_create_histogram():
    pass


def test_histogram():
    pass
