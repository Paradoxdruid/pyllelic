#!/usr/bin/env python3
"""pytest unit tests for __main__.py."""

# Testing
import argparse
import pytest  # noqa  # pylint: disable=unused-import
import sys
from pathlib import Path

# Module to test
import pyllelic.__main__ as main

# mocker.patch('argparse.ArgumentParser.parse_args',
#             return_value=argparse.Namespace(kwarg1=value, kwarg2=value))

EXPECTED_ARGS = argparse.Namespace(
    output_fname="my_data",
    fastq="fh_cellline_tissue.fastq.gz",
    genome="hg19chr5",
    chrom="chr5",
    start=1293000,
    end=1296000,
    viz="plotly",
    fname_pattern="^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$",
    testdir="test",
)
TEST_ARGS = "__main__.py -o my_data -f fh_cellline_tissue.fastq.gz -g hg19chr5 \
                        -chr chr5 -s 1293000 -e 1296000 --viz plotly".split()


def test__parsing(mocker):
    mocker.patch.object(sys, "argv", TEST_ARGS)
    EXPECTED = EXPECTED_ARGS
    actual = main._parsing()

    assert actual == EXPECTED


def test__process_files(mocker):
    mock_process = mocker.patch("pyllelic.__main__.process")

    main._process_files(EXPECTED_ARGS)

    mock_process.prepare_genome.assert_called_with(Path.cwd() / "hg19chr5")
    mock_process.bismark.assert_called_with(
        Path.cwd() / "hg19chr5", Path.cwd() / "fh_cellline_tissue.fastq.gz"
    )
    mock_process.pysam_sort.assert_called_with(
        Path.cwd() / "fh_cellline_tissue.fastq.bam"
    )
    mock_process.pysam_index.assert_called_with(
        Path.cwd() / "fh_cellline_tissue.fastq_sorted.bam"
    )


def test__call_pyllelic(mocker):
    mock_pyllelic = mocker.patch("pyllelic.__main__.pyllelic")

    main._call_pyllelic(EXPECTED_ARGS)

    mock_pyllelic.configure.assert_called_with(
        base_path=str(Path.cwd()),
        prom_file="genome.txt",
        prom_start=1293000,
        prom_end=1296000,
        chrom="chr5",
        offset=1293000,
        viz_backend="plotly",
        fname_pattern=r"^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$",
        test_dir="test",
    )

    mock_pyllelic.make_list_of_bam_files.assert_called_once()

    mock_pyllelic.pyllelic.assert_called_once()


def test_run_pyllelic(mocker, capsys):
    mocker.patch.object(sys, "argv", TEST_ARGS)
    mock_process = mocker.patch("pyllelic.__main__.process")
    mock_pyllelic = mocker.patch("pyllelic.__main__.pyllelic")

    main.run_pyllelic()

    mock_process.prepare_genome.assert_called_with(Path.cwd() / "hg19chr5")
    mock_process.bismark.assert_called_with(
        Path.cwd() / "hg19chr5", Path.cwd() / "fh_cellline_tissue.fastq.gz"
    )
    mock_process.pysam_sort.assert_called_with(
        Path.cwd() / "fh_cellline_tissue.fastq.bam"
    )
    mock_process.pysam_index.assert_called_with(
        Path.cwd() / "fh_cellline_tissue.fastq_sorted.bam"
    )

    mock_pyllelic.configure.assert_called_with(
        base_path=str(Path.cwd()),
        prom_file="genome.txt",
        prom_start=1293000,
        prom_end=1296000,
        chrom="chr5",
        offset=1293000,
        viz_backend="plotly",
        fname_pattern=r"^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$",
        test_dir="test",
    )

    mock_pyllelic.make_list_of_bam_files.assert_called_once()

    mock_pyllelic.pyllelic.assert_called_once()

    # mock_pyllelic.GenomicPositionData.__init__.assert_called_once()

    captured = capsys.readouterr()

    assert (
        captured.out
        == """Preparing genome and processing fastq file... \
this will take a while.\nRunning pyllelic...\nData processed, \
saving output files...\nPyllelic run complete\n"""
    )
