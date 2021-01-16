#!/usr/bin/env python3

import pytest  # noqa
import pandas as pd
import pyllelic


def test_genome_range():
    gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    result = pyllelic.genome_range(position=2, genome_string=gen_str, offset=40)
    assert result == gen_str[8:37]
    assert type(result) == str


def test_samtools_index():
    pass


def test_extract_cell_types():
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
    # input_list_of_dicts = {"TEST11": pd.DataFrame({})}
    pass


def test_return_read_values():
    in_pos = 1
    in_key = "TEST1"
    in_dict1 = {}
    in_dict2 = {}
    in_dict3 = {}
    input_list_of_dicts = {
        "TEST11": pd.DataFrame(in_dict1),
        "TEST2": pd.DataFrame(in_dict2),
        "TEST3": pd.DataFrame(in_dict3),
    }
    in_min_reads = 4
    in_min_sites = 2

    expected = {}

    result = pyllelic.return_read_values(
        pos=in_pos,
        key=in_key,
        dict_of_dfs=input_list_of_dicts,
        min_reads=in_min_reads,
        min_sites=in_min_sites,
    )

    assert result == expected


# def test_simple_pass():
#     """The simplest passing test"""
#     pass


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
