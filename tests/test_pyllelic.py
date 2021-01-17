#!/usr/bin/env python3
"""pytest unit tests for pyllelic."""

import pytest  # noqa
import pandas as pd
import pyllelic
import numpy as np

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


def test_genome_range():
    """Check if correct genome string is returned."""
    gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    result = pyllelic.genome_range(position=2, genome_string=gen_str, offset=40)
    assert result == gen_str[8:37]
    assert isinstance(result, str)


def test_samtools_index():
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


def test_return_individual_data():
    """Check whether the expected and result DataFrames are identical."""
    in_pos = ["1", "2", "3"]
    in_cell = ["TEST1", "TEST2", "TEST3"]
    result = pyllelic.return_individual_data(
        dict_of_dfs=SAMPLE_DICT_OF_DFS, positions=in_pos, cell_types=in_cell
    )

    intermediate = pd.DataFrame.from_dict(
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
            "TEST3": [
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                [1.0, 1.0, 1.0, (2 / 3), (2 / 3), np.nan],
            ],
        },
        orient="index",
        columns=["1", "2", "3"],
    )

    expected = intermediate.astype("object")

    pd.testing.assert_frame_equal(result, expected)


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


def test_write_bam_output_files(tmp_path):
    """Check if bam outputs would be correctly written."""
    pyllelic.config.base_directory = tmp_path
    test_positions = ["1", "2", "3"]
    test_df = pd.DataFrame.from_dict(
        {
            "TEST1": [],
        },
    )
    test_sams = tmp_path / "fh_TEST1_CELL.TERT.bam"

    pyllelic.write_bam_output_files(
        sams=test_sams, positions=test_positions, df=test_df
    )

    expected = tmp_path / "bam_output" / "fh_TEST1_CELL.TERT.bam"

    assert expected.isdir()


# def test_simple_assertions():
#     """Demonstrates passing tests that use assert"""
#     assert True
#     assert [1]
#     assert dict(pytest="awesome")


# def test_negative_assertions():
#     """Demonstrates passing tests that use negated assertions"""
#     assert not False
#     assert not []
#     assert not dict()


# def test_expected_exception():
#     """Demonstrates pytest's raises context manager"""

#     with pytest.raises(ZeroDivisionError):
#         1 / 0

#     with pytest.raises(IOError):
#         open("/some/bogus/file.txt")
