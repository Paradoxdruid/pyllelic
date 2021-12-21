#!/usr/bin/env python3
"""pytest unit tests for quma."""

# Testing
import pytest  # noqa

# Module to test
import pyllelic.quma as quma

# Common test data
# EXPECTED_ALIGN_MATCH = {
#     "qAli": "ATCGATCCGGCATACG",
#     "gAli": "ATCGATCCGGCATACG",
#     "gap": 0,
#     "menum": 3,
#     "unconv": 0,
#     "conv": 3,
#     "pconv": 100.0,
#     "match": 16,
#     "val": "111",
#     "perc": 100.0,
#     "aliMis": 0,
#     "aliLen": 16,
# }

EXPECTED_ALIGN_MATCH = quma.Result(
    qAli="ATCGATCCGGCATACG",
    gAli="ATCGATCCGGCATACG",
    gap=0,
    menum=3,
    unconv=0,
    conv=3,
    pconv=100.0,
    match=16,
    val="111",
    perc=100.0,
    aliMis=0,
    aliLen=16,
)

EXPECTED_ALIGN_MISMATCH = quma.Result(
    qAli="ATTGATCCGGCATACG",
    gAli="ATCGATCCGGCATACG",
    gap=0,
    menum=2,
    unconv=1,
    conv=2,
    pconv=66.7,
    match=17,
    val="011",
    perc=106.2,
    aliMis=-1,
    aliLen=16,
)

TEST_FORMAT_DICT = quma.Result(
    qAli="ATCGATCCGGCATACG",
    gAli="ATCGATCCGGCATACG",
    aliLen=16,
    aliMis=0,
    perc=100.0,
    gap=0,
    menum=0,
    unconv=0,
    conv=3,
    pconv=100.0,
    val=100.0,
)

TEST_SUMMARY_REF = quma.Result(
    qAli="ATCGATCCGGCATACG",
    gAli="ATCGATCCGGCATACG",
    gap=0,
    menum=3,
    unconv=0,
    conv=3,
    pconv=0,
    match=16,
    val="111",
    perc=0,
    aliMis=0,
    aliLen=16,
)

TEST_ALIGN_REF = quma.Result(
    qAli="ATCGATCCGGCATACG",
    gAli="ATCGATCCGGCATACG",
    gap=0,
    menum=0,
    unconv=0,
    conv=0,
    pconv="",
    match=0,
    val="",
    perc="",
    aliMis=0,
)

TEST_ALIGN_REF_MISMATCH = quma.Result(
    qAli="ATTGATCCGGCATACG",
    gAli="ATCGATCCGGCATACG",
    gap=0,
    menum=0,
    unconv=0,
    conv=0,
    pconv="",
    match=0,
    val="",
    perc="",
    aliMis=0,
)


class Test_Quma:
    """Class to test Quma object initialization and read."""

    # pylint: disable=no-self-use

    @pytest.fixture()
    def set_up_quma(self):
        TEST_GSEQ = ">query\nATCGTAGTCGA"
        TEST_QSEQ = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
        return quma.Quma(TEST_GSEQ, TEST_QSEQ)

    def test_init(self, set_up_quma):
        quma_result = set_up_quma
        actual = quma_result.values

        # EXPECTED = (
        #     "genome\t0\tATCGTAGTCGA\t1\t0\n"
        #     + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
        #     + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
        #     + "2\tquery2\tATCGATAGCATT\tATCG-TAGTCGA\tATCGATAGCATT\t"
        #     + "12\t4\t66.7\t1\t1\t0\t1\t100.0\t1\t1\t1\n"
        # )

        EXPECTED = (
            "genome\t0\tATCGTAGTCGA\t1\t0\n"
            + "1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t"
            + "11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n"
            + "2\tquery2\tATCGATAGCATT\tATCG-TAGT\tATCGATAGC\t"
            + "9\t1\t88.9\t1\t1\t0\t1\t100.0\t1\t1\t1\n"
        )

        assert EXPECTED == actual

    def test__parse_genome(self, set_up_quma):
        quma_result = set_up_quma
        EXPECTED = "ATCGTAGTCGA"
        actual = quma_result._parse_genome()
        assert EXPECTED == actual

    def test__parse_biseq(self, set_up_quma):
        quma_result = set_up_quma
        EXPECTED = [
            quma.Fasta(com="query1", pos=None, seq="ATCGTAGTCGA"),
            quma.Fasta(com="query2", pos=None, seq="ATCGATAGCATT"),
        ]
        quma_result._qfile_contents = ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
        actual = quma_result._parse_biseq()
        assert EXPECTED == actual

    def test__parse_biseq_empty_line(self, set_up_quma):
        quma_result = set_up_quma
        EXPECTED = [
            quma.Fasta(com="query1", pos=None, seq=""),
            quma.Fasta(com="query2", pos=None, seq="ATCGATAGCATT"),
        ]
        quma_result._qfile_contents = ">query1\nEEEEEEE\n>query2\nATCGATAGCATT"
        actual = quma_result._parse_biseq()
        assert EXPECTED == actual

    def test__parse_seq(self, set_up_quma):
        quma_result = set_up_quma
        TEST_SEQ = ">query\nATCGTAGTCGA"
        EXPECTED = "ATCGTAGTCGA"
        actual = quma_result._parse_seq(TEST_SEQ)
        assert EXPECTED == actual

    def test__check_char_in_allowed(self, set_up_quma):
        quma_result = set_up_quma
        SEQ = "ABCDEFGHIJ"
        PATTERN = "JABC"

        EXPECTED = "ABCJ"
        actual = quma_result._check_char_in_allowed(SEQ, PATTERN)
        assert EXPECTED == actual

    def test__fasta_make(self, set_up_quma):
        quma_result = set_up_quma
        TEST_SEQ = "ATCGTAGTCGA"
        TEST_SEQ_NAME = "read1"
        EXPECTED = ">read1\nATCGTAGTCGA\n"
        actual = quma_result._fasta_make(TEST_SEQ, TEST_SEQ_NAME)
        assert EXPECTED == actual

    def test__process_fasta_output(self, set_up_quma):
        quma_result = set_up_quma
        TEST_QSEQ = [
            quma.Fasta(com="read0", seq="ATCGATCCGGCATACG"),
            quma.Fasta(com="read1", seq="ATCGATCCGGCATACG"),
        ]
        TEST_QFILEF = "queryF"
        TEST_QFILER = "queryR"
        TEST_GFILEPF = ">genomeF\nATCGATCCGGCATACG"
        # TEST_CPGF = {"2": 1, "7": 1, "14": 1}
        EXPECTED = [
            quma.Reference(
                fasta=quma.Fasta(com="read0", seq="ATCGATCCGGCATACG", pos="1"),
                res=EXPECTED_ALIGN_MATCH,
                dir=1,
                gdir=1,
                exc=1,
            ),
            quma.Reference(
                fasta=quma.Fasta(com="read1", seq="ATCGATCCGGCATACG", pos="2"),
                res=EXPECTED_ALIGN_MATCH,
                dir=1,
                gdir=1,
                exc=1,
            ),
        ]

        actual = quma_result._process_fasta_output(
            TEST_QSEQ,
            TEST_QFILEF,
            TEST_QFILER,
            TEST_GFILEPF,
        )
        assert EXPECTED == actual

    # @pytest.mark.parametrize(
    #     "TEST_SEQ, EXPECTED",
    #     [("ATCGTAGTCGA", "TCGACTACGAT"), ("ATCGTAGTCGO", "OCGACTACGAT")],
    # )
    def test_rev_comp(self, set_up_quma):
        TEST_SEQ = "ATCGTAGTCGA"
        EXPECTED = "TCGACTACGAT"
        quma_result = set_up_quma
        actual = quma_result._rev_comp(TEST_SEQ)
        assert EXPECTED == actual

    def test__align_seq_and_generate_stats(self, set_up_quma):
        quma_result = set_up_quma
        TEST_GFILE = ">genome\nATCGATCCGGCATACG\n"
        TEST_QFILE = ">read1\nATCGATCCGGCATACG\n"
        # TEST_CPG = {"2": 1, "7": 1, "14": 1}
        EXPECTED = EXPECTED_ALIGN_MATCH
        actual = quma_result._align_seq_and_generate_stats(TEST_GFILE, TEST_QFILE)
        assert EXPECTED == actual

    def test__process_alignment_matches(self, set_up_quma):
        quma_result = set_up_quma
        TEST_REF = TEST_ALIGN_REF
        EXPECTED = EXPECTED_ALIGN_MATCH
        actual = quma_result._process_alignment_matches(TEST_REF)
        assert EXPECTED == actual

    def test__process_alignment_matches_mismatch(self, set_up_quma):
        quma_result = set_up_quma
        TEST_REF = TEST_ALIGN_REF_MISMATCH
        EXPECTED = EXPECTED_ALIGN_MISMATCH

        actual = quma_result._process_alignment_matches(TEST_REF)
        assert EXPECTED == actual

    def test__generate_summary_stats(self, set_up_quma):
        quma_result = set_up_quma
        TEST_REF = TEST_SUMMARY_REF
        EXPECTED = EXPECTED_ALIGN_MATCH
        actual = quma_result._generate_summary_stats(TEST_REF)
        assert EXPECTED == actual

    def test__generate_summary_stats_bad(self, set_up_quma):
        quma_result = set_up_quma
        TEST_REF = TEST_SUMMARY_REF
        TEST_REF.conv = 0
        EXPECTED = EXPECTED_ALIGN_MATCH
        EXPECTED.pconv = 0
        EXPECTED.conv = 0
        actual = quma_result._generate_summary_stats(TEST_REF)
        assert EXPECTED == actual

    def test__percentage(self, set_up_quma):
        quma_result = set_up_quma
        TEST_SUM_A = 3
        TEST_SUM_B = 7
        EXPECTED_SUM = 30.0
        actual_sum = quma_result._percentage(TEST_SUM_A, TEST_SUM_B, calc_type="sum")
        assert EXPECTED_SUM == actual_sum

        TEST_TOTAL_A = 3
        TEST_TOTAL_B = 6
        EXPECTED_TOTAL = 50.0
        actual_total = quma_result._percentage(
            TEST_TOTAL_A, TEST_TOTAL_B, calc_type="total"
        )
        assert EXPECTED_TOTAL == actual_total

    def test__percentage_invalid(self, set_up_quma):
        quma_result = set_up_quma
        TEST_SUM_A = 3
        TEST_SUM_B = 7
        with pytest.raises(ValueError):
            _ = quma_result._percentage(TEST_SUM_A, TEST_SUM_B, calc_type="FAKE")

    def test__find_best_dataset(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES.aliMis = 2
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)

        EXPECTED_RES.aliMis = 2

        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_rev(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FFRES.aliMis = 2
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_RES.aliMis = 0
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_fwd_perc(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES.perc = 70.0
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_RES.perc = 70.0
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_rev_perc(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FFRES.perc = 100.0
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_RES.perc = 100.0
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_fwd_unconv(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES.unconv = 1
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_RES.perc = 100.0
        EXPECTED_RES.unconv = 1
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_rev_unconv(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FFRES.unconv = 1
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_RES.unconv = 0
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_fwd_pconv(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES.pconv = 70.0
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_RES.pconv = 70.0
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__find_best_dataset_rev_pconv(self, set_up_quma):
        quma_result = set_up_quma
        TEST_FFRES = EXPECTED_ALIGN_MATCH
        TEST_FRRES = EXPECTED_ALIGN_MATCH
        TEST_FFRES.pconv = 70.0
        EXPECTED_RES = EXPECTED_ALIGN_MATCH
        EXPECTED_DIRECTION = -1
        actual_res, actual_dir = quma_result._find_best_dataset(TEST_FFRES, TEST_FRRES)
        assert EXPECTED_RES == actual_res
        assert EXPECTED_DIRECTION == actual_dir

    def test__format_output(self, set_up_quma):
        quma_result = set_up_quma
        TEST_SEQ = "ATCGATCCGGCATACG"
        TEST_DATA = [
            quma.Reference(
                fasta=quma.Fasta(pos="0", com="read0", seq="ATCGATCCGGCATACG"),
                res=TEST_FORMAT_DICT,
                dir=0,
                gdir=0,
                exc=0,
            ),
            quma.Reference(
                fasta=quma.Fasta(pos="0", com="read1", seq="ATCGATCCGGCATACG"),
                res=TEST_FORMAT_DICT,
                dir=0,
                gdir=0,
                exc=0,
            ),
        ]

        EXPECTED = (
            "genome\t0\tATCGATCCGGCATACG\t1\t0\n0\tread0\t"
            + "ATCGATCCGGCATACG\tATCGATCCGGCATACG\tATCGATCCGGCATACG\t16\t0\t100.0"
            + "\t0\t0\t0\t3\t100.0\t100.0\t0\t0\n0\tread1\tATCGATCCGGCATACG\t"
            + "ATCGATCCGGCATACG\tATCGATCCGGCATACG\t16\t0\t100.0\t0\t0\t0\t3\t"
            + "100.0\t100.0\t0\t0\n"
        )
        actual = quma_result._format_output(TEST_SEQ, TEST_DATA)
        print(repr(actual))
        assert EXPECTED == actual
