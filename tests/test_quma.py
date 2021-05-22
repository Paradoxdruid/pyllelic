#!/usr/bin/env python3
"""pytest unit tests for quma."""

# Testing
import pytest  # noqa
import unittest.mock as mock
import os
import tempfile
from contextlib import contextmanager

# Required libraries for test data

# Module to test
import pyllelic.quma as quma


@contextmanager
def tempinput(data):
    """Helper for virtual files."""
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.write(data)
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)


def test_check_char_in_allowed():
    """Test if string with unallowed characters is returned."""
    SEQ = "ABCDEFGHIJ"
    PATTERN = "JABC"

    EXPECTED = "ABCJ"
    actual = quma.check_char_in_allowed(SEQ, PATTERN)
    assert EXPECTED == actual


def test_curate_seq():
    """Test if string with unallowed characters is returned."""
    SEQ = "ABCDEFGHIJ"

    EXPECTED = "ABCDGH"
    actual = quma.curate_seq(SEQ)
    assert EXPECTED == actual


def test_parse_seq():
    TEST_SEQ = ">query\nATCGTAGTCGA"
    EXPECTED = "ATCGTAGTCGA"
    actual = quma.parse_seq(TEST_SEQ)
    assert EXPECTED == actual


def test_parse_genome():
    TEST_SEQ = b">query\nATCGTAGTCGA"
    EXPECTED = "ATCGTAGTCGA"
    with tempinput(TEST_SEQ) as tempfilename:
        actual = quma.parse_genome(tempfilename)
    assert EXPECTED == actual


def test_multi_fasta_parse():
    TEST_SEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
    EXPECTED = [
        {"com": "query1", "seq": "ATCGTAGTCGA"},
        {"com": "query2", "seq": "ATCGATAGCATT"},
    ]
    actual = quma.multi_fasta_parse(TEST_SEQ)
    assert EXPECTED == actual


def test_parse_biseq():
    TEST_SEQ = b">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
    EXPECTED = [
        {"com": "query1", "seq": "ATCGTAGTCGA"},
        {"com": "query2", "seq": "ATCGATAGCATT"},
    ]
    with tempinput(TEST_SEQ) as tempfilename:
        actual = quma.parse_biseq(tempfilename)
    assert EXPECTED == actual


def test_parse_multi():
    TEST_SEQ = b">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
    EXPECTED = (
        None,
        [
            {"com": "query1", "seq": "ATCGTAGTCGA"},
            {"com": "query2", "seq": "ATCGATAGCATT"},
        ],
    )
    with tempinput(TEST_SEQ) as tempfilename:
        actual = quma.parse_multi(tempfilename)
    assert EXPECTED == actual


def test_fasta_make():
    TEST_SEQ = "ATCGTAGTCGA"
    TEST_SEQ_NAME = "read1"
    EXPECTED = ">read1\nATCGTAGTCGA\n"
    actual = quma.fasta_make(TEST_SEQ, TEST_SEQ_NAME)
    assert EXPECTED == actual


def test_fasta_print(mocker):
    open_mock = mocker.mock_open()
    TEST_SEQ = "ATCGTAGTCGA"
    TEST_SEQ_NAME = "read1"
    TEST_PATH = "/Users/user/test"
    with mock.patch("builtins.open", open_mock, create=True) as m:
        quma.fasta_print(TEST_SEQ, TEST_SEQ_NAME, TEST_PATH)

    open_mock.assert_called_with(TEST_PATH, "w")
    handle = m()
    handle.write.assert_called_with(">read1\nATCGTAGTCGA\n")


def test_fasta_output():
    TEST_SEQ = "ATCGTAGTCGA"
    TEST_SEQ_NAME = "read1"
    EXPECTED = ">read1\nATCGTAGTCGA\n"
    actual = quma.fasta_output(TEST_SEQ, TEST_SEQ_NAME)
    assert EXPECTED == actual


def test_rev_comp():
    pass


def test_align_seq_and_generate_stats():
    pass


def test__generate_summary_stats():
    pass


def test__percentage():
    pass


def test_process_alignment_matches():
    pass


def test_process_fasta_output():
    pass


def test_format_output():
    pass


def test_find_cpg():
    pass


def test_quma_main():
    pass
