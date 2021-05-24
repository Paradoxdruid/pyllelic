#!/usr/bin/env python3
"""pytest unit tests for process."""

# Testing
import pytest  # noqa  # pylint: disable=unused-import
import unittest.mock as mock
import os
import tempfile
from contextlib import contextmanager

# Required libraries for test data
from pathlib import Path

# Module to test
import pyllelic.process as process  # noqa  # pylint: disable=unused-import


@contextmanager
def tempinput(data, suffix):  # pragma: no cover
    """Helper for virtual files."""
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp.write(data)
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)


def test_process_fastq_to_list():
    pass


def test_make_records_to_dictionary():
    pass


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
        TEST_FASTQ.stem + ".bam",
    ]
    _ = process.bowtie2_fastq_to_bam(TEST_INDEX, TEST_FASTQ, TEST_CORES)

    mock_subp.run.assert_called_once_with(
        TEST_COMMAND, capture_output=True, text=True, check=True
    )


@mock.patch("pyllelic.process.pysam")
def test_process_pysam_sort(mock_pysam):
    TEST_PATH = Path().cwd()
    _ = process.process_pysam_sort(TEST_PATH)
    mock_pysam.sort.assert_called_once_with(
        "-o", f">{TEST_PATH.stem}_sorted.bam", os.fspath(TEST_PATH)
    )


@mock.patch("pyllelic.process.pysam")
def test_pysam_index(mock_pysam):
    TEST_PATH = Path().cwd()
    _ = process.process_pysam_index(TEST_PATH)
    mock_pysam.index.assert_called_once_with(os.fspath(TEST_PATH))
