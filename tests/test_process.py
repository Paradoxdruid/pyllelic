#!/usr/bin/env python3
"""pytest unit tests for process."""

# Testing
import pytest  # noqa  # pylint: disable=unused-import
import unittest.mock as mock
import os
import tempfile

# Required libraries for test data
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

# Module to test
import pyllelic.process as process  # noqa  # pylint: disable=unused-import


# Constants
EXPECTED_SEQ_RECORD = SeqIO.SeqRecord(
    seq=Seq("GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT"),
    id="SEQ_ID",
    name="SEQ_ID",
    description="SEQ_ID",
    dbxrefs=[],
)

TEST_FASTQ_CONTENTS = b"""@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65"""


def test_process_fastq_to_list():
    EXPECTED = EXPECTED_SEQ_RECORD
    with tempfile.NamedTemporaryFile(suffix=".fastq", prefix="test1") as my_file:
        FASTQ_CONTENTS = TEST_FASTQ_CONTENTS
        my_file.write(FASTQ_CONTENTS)
        my_file.seek(0)
        TEST_FILEPATH = Path(my_file.name)
        actual = process.fastq_to_list(TEST_FILEPATH)

    assert EXPECTED.seq == actual[0].seq


def test_process_fastq_to_list_wrong_filetype():
    with tempfile.NamedTemporaryFile(suffix=".txt", prefix="test1") as my_file:
        TEST_FILEPATH = Path(my_file.name)
        actual = process.fastq_to_list(TEST_FILEPATH)

    assert actual is None


def test_process_fastq_to_list_wrong_filetype_multi_extension():
    with tempfile.NamedTemporaryFile(suffix=".fastq.txt", prefix="test1") as my_file:
        TEST_FILEPATH = Path(my_file.name)
        actual = process.fastq_to_list(TEST_FILEPATH)

    assert actual is None


def test_process_fastq_to_list_gz():
    EXPECTED = EXPECTED_SEQ_RECORD
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", prefix="test1") as my_file:
        FASTQ_CONTENTS = TEST_FASTQ_CONTENTS
        gzip_file = gzip.GzipFile(mode="wb", fileobj=my_file)
        gzip_file.write(FASTQ_CONTENTS)
        gzip_file.close()
        my_file.seek(0)
        TEST_FILEPATH = Path(my_file.name)
        actual = process.fastq_to_list(TEST_FILEPATH)

    assert EXPECTED.seq == actual[0].seq


def test_make_records_to_dictionary():
    TEST_RECORD_LIST = [
        SeqIO.SeqRecord(
            Seq("ATGCTCGTAGCTGATCGA"),
            id="test1",
            name="test1",
            description="test record #1",
        ),
        SeqIO.SeqRecord(
            Seq("GTGCTCGTAGCTGATCGA"),
            id="test2",
            name="test2",
            description="test record #2",
        ),
    ]
    EXPECTED = {"test1": TEST_RECORD_LIST[0], "test2": TEST_RECORD_LIST[1]}
    actual = process.make_records_to_dictionary(TEST_RECORD_LIST)
    assert EXPECTED == actual


@mock.patch("pyllelic.process.subprocess")
def test_build_bowtie2_index(mock_subp):
    TEST_FASTQ = Path("/Users/user/test.fastq")
    TEST_COMMAND = ["bowtie2-build", "index", os.fspath(TEST_FASTQ)]

    _ = process.build_bowtie2_index(TEST_FASTQ)

    mock_subp.run.assert_called_once_with(
        TEST_COMMAND, capture_output=True, text=True, check=True
    )


@mock.patch("pyllelic.process.subprocess")
def test_bowtie2_fastq_to_bam(mock_subp):
    TEST_CORES = 4
    TEST_INDEX = Path("/Users/user/bowtie_index")
    TEST_FASTQ = Path("/Users/user/test.fastq")
    TEST_COMMAND = [
        "bowtie2",
        "-p",
        str(TEST_CORES),
        "-x",
        os.fspath(TEST_INDEX),
        "-U",
        os.fspath(TEST_FASTQ),
        "|",
        "samtools",
        "view",
        "-bS",
        "-",
        ">",
        "/Users/user/" + TEST_FASTQ.stem + ".bam",
    ]
    _ = process.bowtie2_fastq_to_bam(TEST_INDEX, TEST_FASTQ, TEST_CORES)

    mock_subp.run.assert_called_once_with(
        TEST_COMMAND, capture_output=True, text=True, check=True
    )


@mock.patch("pyllelic.process.pysam")
def test_process_pysam_sort(mock_pysam):
    TEST_PATH = Path().cwd()
    _ = process.pysam_sort(TEST_PATH)
    mock_pysam.sort.assert_called_once_with(
        "-o", f"{TEST_PATH.parent}/{TEST_PATH.stem}_sorted.bam", os.fspath(TEST_PATH)
    )


@mock.patch("pyllelic.process.pysam")
def test_pysam_index(mock_pysam):
    TEST_PATH = Path().cwd()
    _ = process.pysam_index(TEST_PATH)
    mock_pysam.index.assert_called_once_with(os.fspath(TEST_PATH))


@mock.patch("pyllelic.process.subprocess")
def test_prepare_genome(mock_subp):
    TEST_INDEX = Path("/Users/user/bowtie_index")
    TEST_ALIGNER = Path("/usr/bin/bowtie2/")
    TEST_COMMAND = [
        "bismark_genome_preparation",
        "--path_to_aligner",
        str(TEST_ALIGNER),
        str(TEST_INDEX),
    ]
    _ = process.prepare_genome(TEST_INDEX, TEST_ALIGNER)

    mock_subp.run.assert_called_once_with(
        TEST_COMMAND,
        capture_output=True,
        text=True,
        check=True,
        cwd=TEST_INDEX.parent,
    )


@mock.patch("pyllelic.process.subprocess")
def test_prepare_genome_no_aligner(mock_subp):
    TEST_INDEX = Path("/Users/user/bowtie_index")
    TEST_COMMAND = [
        "bismark_genome_preparation",
        str(TEST_INDEX),
    ]
    _ = process.prepare_genome(TEST_INDEX)

    mock_subp.run.assert_called_once_with(
        TEST_COMMAND,
        capture_output=True,
        text=True,
        check=True,
        cwd=TEST_INDEX.parent,
    )


@mock.patch("pyllelic.process.subprocess")
def test_bismark(mock_subp):
    TEST_GENOME = Path("/Users/user/bowtie_index")
    TEST_FASTQ = Path("/Users/user/test.fastq")
    TEST_COMMAND = [
        "bismark",
        "--genome",
        str(TEST_GENOME),
        str(TEST_FASTQ),
    ]
    _ = process.bismark(TEST_GENOME, TEST_FASTQ)

    mock_subp.run.assert_called_once_with(
        TEST_COMMAND,
        capture_output=True,
        text=True,
        check=True,
        cwd=TEST_FASTQ.parent,
    )


def test_retrieve_promoter_seq(requests_mock):
    EXPECTED = (
        "cgcgtgtccatcaaaacgtgaaggtgaacctcgtaagtttatgcaaactggacaggagggagagcaga"
        + "ggcagagatcaccgtgtccactcgacgtcctgagcgaaaagccacgtgtgcccacgtgacgatggagac"
        + "aggaggaccagggctctgcctgcccccttttctgagcccctactgcattcagctctggggcctg"
    )
    requests_mock.get(
        "https://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr5:1293200,1293400",
        text='<?xml version="1.0" standalone="no"?>\n<!DOCTYPE DASDNA SYSTEM '
        + '"http://www.biodas.org/dtd/dasdna.dtd">\n<DASDNA>\n<SEQUENCE id="chr5" '
        + 'start="1293200" stop="1293400" version="1.00">\n<DNA length="201">\n'
        + "cgcgtgtccatcaaaacgtgaaggtgaacctcgtaagtttatgcaaactg\ngacaggagggagagcaga"
        + "ggcagagatcaccgtgtccactcgacgtcctg\nagcgaaaagccacgtgtgcccacgtgacgatggagac"
        + "aggaggaccaggg\nctctgcctgcccccttttctgagcccctactgcattcagctctggggcct\ng\n"
        + "</DNA>\n</SEQUENCE>\n</DASDNA>\n",
    )

    with tempfile.NamedTemporaryFile(suffix=".txt", prefix="promoter") as my_file:
        process.retrieve_promoter_seq(my_file.name, "chr5", 1293200, 1293400)
        actual = Path(my_file.name).read_text()
        assert actual == EXPECTED
