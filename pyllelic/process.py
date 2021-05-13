#!/usr/bin/env python3
"""Utilities to pre-process and prepare data for use in pyllelic.
"""

import gzip
from Bio import SeqIO
from pathlib import Path
from typing import List, Dict, Optional
import subprocess
import os


def process_fastq_to_list(filepath: Path) -> Optional[List[SeqIO.SeqRecord]]:
    """Read a .fastq or fastq.gz file into an in-memory record_list.

    This is a time and memory intensive operation!

    Args:
        filepath (Path): file path to a fastq.gz file

    Returns:
        List[SeqRecord]: list of biopython sequence records from the fastq file
    """

    if ".fastq" not in filepath.suffixes:
        print("Wrong filetype")
        return None

    record_list: List[SeqIO.SeqRecord] = []
    if ".gz" in filepath.suffixes[-1]:
        with gzip.open(filepath, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                record_list.append(record)
        return record_list

    if ".fastq" in filepath.suffixes[-1]:
        with open(filepath, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                record_list.append(record)
        return record_list

    # If doesn't match readable suffixes
    print("Wrong filetype")
    return None


def make_records_to_dictionary(
    record_list: List[SeqIO.SeqRecord],
) -> Dict[str, SeqIO.SeqRecord]:
    """Take in list of biopython SeqRecords and output a dictionary
       with keys of the record name.

    Args:
        record_list (List[SeqRecord]): biopython sequence records from a fastq file

    Returns:
        Dict[str, SeqRecord]: dict of biopython SeqRecords from a fastq file
    """
    return dict(zip([record.id for record in record_list], record_list))


def build_bowtie2_index(fasta: Path) -> str:
    """Helper function to run external bowtie2-build tool.

    Args:
        fasta (Path): filepath to fasta file to build index from

    Returns:
        str: output from bowtie2-build shell command, usually discarded
    """

    command: List[str] = ["bowtie2-build", "index", os.fspath(fasta)]

    output: subprocess.CompletedProcess = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: str = output.stdout

    return out


def bowtie2_fastq_to_bam(index: Path, fastq: Path, cores: int) -> bytes:
    """Helper function to run external bowtie2-build tool.

    Args:
        index (Path): filepath to bowtie index file
        fastq (Path): filepath to fastq file to convert to bam
        cores (int): number of cores to use for processing

    Returns:
        bytes: output from bowtie2 and samtools shell command, usually discarded
    """

    command: List[str] = [
        "bowtie2",
        "-p",
        str(cores),
        "-x",
        os.fspath(index),
        "-U",
        os.fspath(fastq),
        "|",
        "samtools",
        "view",
        "-bS",
        "-",
        ">",
        fastq.stem + ".bam",
    ]

    output: subprocess.CompletedProcess = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: bytes = output.stdout

    return out


def samtools_sort(bamfile: Path) -> str:
    """Helper function to run external samtools sort.

    Args:
        bamfile (Path): filepath to bam file

    Returns:
        str: output from samtools shell command, usually discarded
    """

    command: List[str] = [
        "samtools",
        "sort",
        os.fspath(bamfile),
        ">",
        bamfile.stem + "_sorted.bam",
    ]

    output: subprocess.CompletedProcess = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: str = output.stdout

    return out


def process_samtools_index(bamfile: Path) -> str:
    """Helper function to run external samtools index.

    Args:
        bamfile (Path): filepath to bam file

    Returns:
        str: output from samtools shell command, usually discarded
    """

    command: List[str] = [
        "samtools",
        "index",
        os.fspath(bamfile),
    ]

    output: subprocess.CompletedProcess = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: str = output.stdout

    return out
