#!/usr/bin/env python3
"""pytest unit tests for process."""

# Testing
import pytest  # noqa  # pylint: disable=unused-import

# import unittest.mock as mock
import os
import tempfile
from contextlib import contextmanager

# Required libraries for test data

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


def test_build_bowtie2_index():
    pass


def test_bowtie2_fastq_to_bam():
    pass


def test_process_pysam_sort():  # not tested
    pass


def test_process_pysam_index():  # not tested
    pass
