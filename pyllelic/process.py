#!/usr/bin/env python3
"""Utilities to pre-process and prepare data for use in pyllelic."""

import gzip
import os
import re
import requests
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import pysam
from Bio import SeqIO


def fastq_to_list(filepath: Path) -> Optional[List[SeqIO.SeqRecord]]:
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
        fastq.parent + fastq.stem + ".bam",
    ]

    output: subprocess.CompletedProcess = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: bytes = output.stdout

    return out


def pysam_sort(bamfile: Path) -> bool:
    """Helper function to run pysam samtools sort.

    Args:
        bamfile (Path): filepath to bam file

    Returns:
        bool: verification of samtools command, usually discarded
    """

    pysam.sort("-o", f"{bamfile.parent}/{bamfile.stem}_sorted.bam", os.fspath(bamfile))

    return True


def pysam_index(bamfile: Path) -> bool:
    """Helper function to run external samtools index.

    Args:
        bamfile (Path): filepath to bam file

    Returns:
        bool: verification of samtools command, usually discarded
    """

    pysam.index(os.fspath(bamfile))

    return True


def retrieve_promoter_seq(filename: str, chrom: str, start: int, end: int) -> None:
    """Retrieve the genomic sequence of interest from UCSC Genome Browser.

    Args:
        filename (Path): path to store genomic sequence
        chrom (str): chromosome of interest, e.g. "chr5"
        start (int): start position for region of interest
        end (int): end position for region of interest
    """
    response = requests.get(
        f"https://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment={chrom}:{start},{end}"
    )
    pattern = re.compile(r"<DNA.*>(.*)<\/DNA>", flags=re.DOTALL)
    match = pattern.findall(response.text)
    seq = match[0].replace("\n", "")

    Path(filename).write_text(seq)
