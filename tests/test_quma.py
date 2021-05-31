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


# Common test data
EXPECTED_ALIGN_MATCH = {
    "qAli": "ATCGATCCGGCATACG",
    "gAli": "ATCGATCCGGCATACG",
    "gap": 0,
    "menum": 3,
    "unconv": 0,
    "conv": 3,
    "pconv": 100.0,
    "match": 16,
    "val": "111",
    "perc": 100.0,
    "aliMis": 0,
    "aliLen": 16,
}

TEST_FORMAT_DICT = {
    "qAli": "ATCGATCCGGCATACG",
    "gAli": "ATCGATCCGGCATACG",
    "aliLen": 16,
    "aliMis": 0,
    "perc": 100.0,
    "gap": 0,
    "menum": 0,
    "unconv": 0,
    "conv": 3,
    "pconv": 100.0,
    "val": 100.0,
}

TEST_SUMMARY_REF = {
    "qAli": "ATCGATCCGGCATACG",
    "gAli": "ATCGATCCGGCATACG",
    "gap": 0,
    "menum": 3,
    "unconv": 0,
    "conv": 3,
    "pconv": 0,
    "match": 16,
    "val": "111",
    "perc": 0,
    "aliMis": 0,
    "aliLen": 16,
}

TEST_ALIGN_REF = {
    "qAli": "ATCGATCCGGCATACG",
    "gAli": "ATCGATCCGGCATACG",
    "gap": 0,
    "menum": 0,
    "unconv": 0,
    "conv": 0,
    "pconv": "",
    "match": 0,
    "val": "",
    "perc": "",
    "aliMis": 0,
}


# Test functions
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


@pytest.mark.parametrize(
    "TEST_SEQ, EXPECTED",
    [("ATCGTAGTCGA", "TCGACTACGAT"), ("ATCGTAGTCGO", "OCGACTACGAT")],
)
def test_rev_comp(TEST_SEQ, EXPECTED):
    actual = quma.rev_comp(TEST_SEQ)
    assert EXPECTED == actual


def test_align_seq_and_generate_stats():
    TEST_GFILE = ">genome\nATCGATCCGGCATACG\n"
    TEST_QFILE = ">read1\nATCGATCCGGCATACG\n"
    TEST_CPG = {"2": 1, "7": 1, "14": 1}
    EXPECTED = EXPECTED_ALIGN_MATCH
    actual = quma.align_seq_and_generate_stats(TEST_GFILE, TEST_QFILE, TEST_CPG)
    assert EXPECTED == actual


def test__generate_summary_stats():
    TEST_REF = TEST_SUMMARY_REF
    EXPECTED = EXPECTED_ALIGN_MATCH
    actual = quma._generate_summary_stats(TEST_REF)
    assert EXPECTED == actual


def test__generate_summary_stats_bad():
    TEST_REF = TEST_SUMMARY_REF.copy()
    TEST_REF["conv"] = 0
    EXPECTED = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED["pconv"] = 0
    EXPECTED["conv"] = 0
    actual = quma._generate_summary_stats(TEST_REF)
    assert EXPECTED == actual


def test__percentage():
    TEST_SUM_A = 3
    TEST_SUM_B = 7
    EXPECTED_SUM = "30.0"
    actual_sum = quma._percentage(TEST_SUM_A, TEST_SUM_B, calc_type="sum")
    assert EXPECTED_SUM == actual_sum

    TEST_TOTAL_A = 3
    TEST_TOTAL_B = 6
    EXPECTED_TOTAL = "50.0"
    actual_total = quma._percentage(TEST_TOTAL_A, TEST_TOTAL_B, calc_type="total")
    assert EXPECTED_TOTAL == actual_total


def test_process_alignment_matches():
    TEST_REF = TEST_ALIGN_REF.copy()
    TEST_CPG = {"2": 1, "7": 1, "14": 1}
    EXPECTED = EXPECTED_ALIGN_MATCH
    actual = quma.process_alignment_matches(TEST_REF, TEST_CPG)
    assert EXPECTED == actual


def test_process_alignment_matches_g_gap():
    TEST_REF = TEST_ALIGN_REF.copy()
    TEST_REF["gAli"] = "-TCGATCCGGCATACG"
    TEST_CPG = {"2": 1, "7": 1, "14": 1}
    EXPECTED = {
        "qAli": "ATCGATCCGGCATACG",
        "gAli": "-TCGATCCGGCATACG",
        "gap": 1,
        "menum": 0,
        "unconv": 0,
        "conv": 0,
        "pconv": 0.0,
        "match": 15,
        "val": "-",
        "perc": 93.8,
        "aliMis": 1,
        "aliLen": 16,
    }
    actual = quma.process_alignment_matches(TEST_REF, TEST_CPG)
    assert EXPECTED == actual


def test_process_alignment_matches_q_gap():
    TEST_REF = TEST_ALIGN_REF.copy()
    TEST_REF["qAli"] = "AT-GATCCGGCATACG"
    TEST_CPG = {"2": 1, "7": 1, "14": 1}
    EXPECTED = {
        "qAli": "AT-GATCCGGCATACG",
        "gAli": "ATCGATCCGGCATACG",
        "gap": 1,
        "menum": 2,
        "unconv": 0,
        "conv": 2,
        "pconv": 100.0,
        "match": 15,
        "val": "-11",
        "perc": 93.8,
        "aliMis": 1,
        "aliLen": 16,
    }
    actual = quma.process_alignment_matches(TEST_REF, TEST_CPG)
    assert EXPECTED == actual


def test_process_alignment_matches_q_is_T():
    TEST_REF = TEST_ALIGN_REF.copy()
    TEST_REF["qAli"] = "ATTGATCCGGCATACG"
    TEST_CPG = {"2": 1, "7": 1, "14": 1}
    EXPECTED = {
        "qAli": "ATTGATCCGGCATACG",
        "gAli": "ATCGATCCGGCATACG",
        "gap": 0,
        "menum": 2,
        "unconv": 1,
        "conv": 2,
        "pconv": 66.7,
        "match": 16,
        "val": "011",
        "perc": 100.0,
        "aliMis": 0,
        "aliLen": 16,
    }
    actual = quma.process_alignment_matches(TEST_REF, TEST_CPG)
    assert EXPECTED == actual


def test_process_alignment_matches_q_is_A():
    TEST_REF = TEST_ALIGN_REF.copy()
    TEST_REF["qAli"] = "ATAGATCCGGCATACG"
    TEST_CPG = {"2": 1, "7": 1, "14": 1}
    EXPECTED = {
        "qAli": "ATAGATCCGGCATACG",
        "gAli": "ATCGATCCGGCATACG",
        "gap": 0,
        "menum": 2,
        "unconv": 0,
        "conv": 2,
        "pconv": 100.0,
        "match": 15,
        "val": "A11",
        "perc": 93.8,
        "aliMis": 1,
        "aliLen": 16,
    }
    actual = quma.process_alignment_matches(TEST_REF, TEST_CPG)
    assert EXPECTED == actual


def test__find_best_dataset():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES["aliMis"] = 2
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = 1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_rev():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FFRES["aliMis"] = 2
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_RES["aliMis"] = 0
    EXPECTED_DIRECTION = -1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_fwd_perc():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES["perc"] = 70.0
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = 1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_rev_perc():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FFRES["perc"] = 70.0
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = -1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_fwd_unconv():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES["unconv"] = 1
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = 1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_rev_unconv():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FFRES["unconv"] = 1
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = -1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_fwd_pconv():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES["pconv"] = 70.0
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = 1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test__find_best_dataset_rev_pconv():
    TEST_FFRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FRRES = EXPECTED_ALIGN_MATCH.copy()
    TEST_FFRES["pconv"] = 70.0
    EXPECTED_RES = EXPECTED_ALIGN_MATCH.copy()
    EXPECTED_DIRECTION = -1
    actual_res, actual_dir = quma._find_best_dataset(TEST_FFRES, TEST_FRRES)
    assert EXPECTED_RES == actual_res
    assert EXPECTED_DIRECTION == actual_dir


def test_process_fasta_output():
    TEST_QSEQ = [
        {"com": "read0", "seq": "ATCGATCCGGCATACG"},
        {"com": "read1", "seq": "ATCGATCCGGCATACG"},
    ]
    TEST_QFILEF = "queryF"
    TEST_QFILER = "queryR"
    TEST_GFILEPF = ">genomeF\nATCGATCCGGCATACG"
    TEST_GFILEPR = ">genomeR\nCGTATGCCGGATCGAT"
    TEST_CPGF = {"2": 1, "7": 1, "14": 1}
    TEST_CPGR = {"12": 1, "7": 1, "0": 1}
    EXPECTED = [
        {
            "fa": {"com": "read0", "seq": "ATCGATCCGGCATACG", "pos": "1"},
            "res": EXPECTED_ALIGN_MATCH,
            "dir": 1,
            "gdir": 1,
            "exc": 1,
        },
        {
            "fa": {"com": "read1", "seq": "ATCGATCCGGCATACG", "pos": "2"},
            "res": EXPECTED_ALIGN_MATCH,
            "dir": 1,
            "gdir": 1,
            "exc": 1,
        },
    ]

    actual = quma.process_fasta_output(
        TEST_QSEQ,
        TEST_QFILEF,
        TEST_QFILER,
        TEST_GFILEPF,
        TEST_GFILEPR,
        TEST_CPGF,
        TEST_CPGR,
    )
    assert EXPECTED == actual


def test_format_output():
    TEST_SEQ = "ATCGATCCGGCATACG"
    TEST_POSITIONS = ["2", "7", "14"]
    TEST_DATA = [
        {
            "fa": {"pos": "0", "com": "read0", "seq": "ATCGATCCGGCATACG"},
            "res": TEST_FORMAT_DICT,
            "dir": "0",
            "gdir": "0",
        },
        {
            "fa": {"pos": "0", "com": "read1", "seq": "ATCGATCCGGCATACG"},
            "res": TEST_FORMAT_DICT,
            "dir": "0",
            "gdir": "0",
        },
    ]
    EXPECTED = (
        "genome\t0\tATCGATCCGGCATACG\t3\t2,7,14\n"
        + "0\tread0\tATCGATCCGGCATACG\tATCGATCCGGCATACG\tATCGATCCGGCATACG\t"
        + "16\t0\t100.0\t0\t0\t0\t3\t100.0\t100.0\t0\t0\n"
        + "0\tread1\tATCGATCCGGCATACG\tATCGATCCGGCATACG\tATCGATCCGGCATACG\t"
        + "16\t0\t100.0\t0\t0\t0\t3\t100.0\t100.0\t0\t0\n"
    )
    actual = quma.format_output(TEST_SEQ, TEST_POSITIONS, TEST_DATA)
    assert EXPECTED == actual


def test_find_cpg():
    TEST_SEQ = "ATCGATCCGGCATACG"
    EXPECTED = (
        ["2", "7", "14"],
        {"2": 1, "7": 1, "14": 1},
        {"12": 1, "7": 1, "0": 1},
    )
    actual = quma.find_cpg(TEST_SEQ)
    assert EXPECTED == actual


def test_quma_main():
    TEST_GSEQ = b">query\nATCGTAGTCGA"
    TEST_QSEQ = b">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
    with tempinput(TEST_GSEQ) as temp_gseq:
        with tempinput(TEST_QSEQ) as temp_qseq:
            actual = quma.quma_main(temp_gseq, temp_qseq)
    EXPECTED = (
        "genome\t0\tATCGTAGTCGA\t2\t2,8\n"
        + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
        + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
        + "2\tquery2\tATCGATAGCATT\tATCG-TAGT\tATCGATAGC\t"
        + "9\t1\t88.9\t1\t1\t0\t1\t100.0\t1\t1\t1\n"
    )
    print(actual)
    assert EXPECTED == actual
