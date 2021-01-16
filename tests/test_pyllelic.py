#!/usr/bin/env python3

import pytest  # noqa

import pyllelic


def test_genome_range():
    gen_str = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    result = pyllelic.genome_range(position=2, genome_string=gen_str, offset=40)
    assert result == gen_str[8:37]
    assert type(result) == str


def test_samtools_index():
    pass


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
