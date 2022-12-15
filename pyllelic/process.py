#!/usr/bin/env python3
"""Utilities to pre-process and prepare data for use in pyllelic."""

import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import pysam
import requests
from Bio import SeqIO

from pyllelic.config import Config
from pyllelic.pyllelic import GenomicPositionData


class ShellCommandError(Exception):
    """Error for shell utilities that aren't installed."""


class FileNameError(Exception):
    """Error for invalid filetypes."""


def fastq_to_list(filepath: Path) -> List[SeqIO.SeqRecord]:
    """Read a .fastq or fastq.gz file into an in-memory record_list.

    This is a time and memory intensive operation!

    Args:
        filepath (Path): file path to a fastq.gz file

    Returns:
        List[SeqRecord]: list of biopython sequence records from the fastq file
    """

    if ".fastq" not in filepath.suffixes:
        raise FileNameError("Wrong filetype")

    record_list: List[SeqIO.SeqRecord] = []
    if ".gz" in filepath.suffixes[-1]:
        with gzip.open(filepath, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                record_list.append(record)
        return record_list

    if ".fastq" in filepath.suffixes[-1]:
        with open(filepath) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                record_list.append(record)
        return record_list

    # If doesn't match readable suffixes
    raise FileNameError("Wrong filetype")


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

    Raises:
        ShellCommandError: bowtie2-build is not installed.
    """

    command: List[str] = ["bowtie2-build", "index", os.fspath(fasta)]

    if shutil.which(command[0]) is None:
        raise ShellCommandError("bowtie2-build is not installed.")

    output: subprocess.CompletedProcess[str] = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: str = output.stdout

    return out


def bowtie2_fastq_to_bam(index: Path, fastq: Path, cores: int) -> str:
    """Helper function to run external bowtie2-build tool.

    Args:
        index (Path): filepath to bowtie index file
        fastq (Path): filepath to fastq file to convert to bam
        cores (int): number of cores to use for processing

    Returns:
        str: output from bowtie2 and samtools shell command, usually discarded

    Raises:
        ShellCommandError: bowtie2 is not installed.
    """

    command: List[str] = [
        "bowtie2",
        "-p",
        str(cores),
        "-x",
        str(index),
        "-U",
        str(fastq),
        "|",
        "samtools",
        "view",
        "-bS",
        "-",
        ">",
        str(fastq.parent) + "/" + str(fastq.stem) + ".bam",
    ]

    if shutil.which(command[0]) is None:
        raise ShellCommandError("bowtie2 is not installed.")

    output: subprocess.CompletedProcess[str] = subprocess.run(
        command, capture_output=True, text=True, check=True
    )
    out: str = output.stdout

    return out


def sort_bam(bamfile: Path) -> bool:
    """Helper function to run pysam samtools sort.

    Args:
        bamfile (Path): filepath to bam file

    Returns:
        bool: verification of samtools command, usually discarded
    """

    pysam.sort(  # type:ignore[attr-defined]
        "-o", f"{bamfile.parent}/{bamfile.stem}_sorted.bam", os.fspath(bamfile)
    )

    return True


def index_bam(bamfile: Path) -> bool:
    """Helper function to run external samtools index.

    Args:
        bamfile (Path): filepath to bam file

    Returns:
        bool: verification of samtools command, usually discarded
    """

    pysam.index(os.fspath(bamfile))  # type:ignore[attr-defined]

    return True


def retrieve_seq(filename: str, chrom: str, start: int, end: int) -> None:
    """Retrieve the genomic sequence of interest from UCSC Genome Browser.

    Args:
        filename (str): path to store genomic sequence
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


def prepare_genome(index: Path, aligner: Optional[Path] = None) -> str:
    """Helper function to run external bismark genome preparation tool.

    Uses genomes from, e.g.: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

    Bismark documentation at:
    https://github.com/FelixKrueger/Bismark/tree/master/Docs

    Args:
        index (Path): filepath to unprocessed genome file.
        aligner (Optional[Path]): filepath to bowtie2 alignment program.

    Returns:
        str: output from genome preparation shell command, usually discarded

    Raises:
        ShellCommandError: bismark_genome_preparation is not installed.
    """
    command: List[str]
    if aligner:
        command = [
            "bismark_genome_preparation",
            "--path_to_aligner",
            str(aligner),
            str(index),
        ]
    else:
        command = [
            "bismark_genome_preparation",
            str(index),
        ]

    if shutil.which(command[0]) is None:
        raise ShellCommandError("bismark_genome_preparation is not installed.")

    output: subprocess.CompletedProcess[str] = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True,
        cwd=index.parent,
    )
    out: str = output.stdout

    return out


def bismark(genome: Path, fastq: Path) -> str:
    """Helper function to run external bismark tool.

    Bismark documentation at:
    https://github.com/FelixKrueger/Bismark/tree/master/Docs

    Args:
        genome (Path): filepath to directory of bismark processed genome files.
        fastq (Path): filepath to fastq file to process.

    Returns:
        str: output from bismark shell command, usually discarded

    Raises:
        ShellCommandError: bismark is not installed.
    """

    command: List[str] = [
        "bismark",
        "--genome",
        str(genome),
        str(fastq),
    ]

    if shutil.which(command[0]) is None:
        raise ShellCommandError("bismark is not installed.")

    output: subprocess.CompletedProcess[str] = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True,
        cwd=fastq.parent,
    )
    out: str = output.stdout

    return out


def methbank_bed_to_indiv_data(
    path: Path, chrom: str, start: int, stop: int, viz: str = "plotly"
) -> GenomicPositionData:
    """Helper function to convert MethBank BED file into individual_data df.

    MethBank: https://ngdc.cncb.ac.cn/methbank/

    Args:
        path (Path): path to MethBank formatted BED file

    Returns:
        GenomicPositionData: mostly complete pyllelic object with data from BED.
    """

    df = pd.read_csv(path, sep="\t")
    df = df[df["#chr"] == chrom]
    df = df[(df["start"] > start) & (df["start"] < stop)]
    df["un_num"] = df["total_num"] - df["methy_num"]
    df["mean"] = df["percent_num"]
    df["mode"] = 0
    df.loc[df["methy_num"] >= df["un_num"], "mode"] = 1.0
    df["diff"] = abs(df["mode"] - df["mean"])
    my_ind_dict = {}
    for _, each in df.iterrows():
        val = []
        val.extend([1.0] * each["methy_num"])
        val.extend([0.0] * each["un_num"])
        my_ind_dict[each.start] = val
    ind_df = pd.DataFrame({str(k): [v] for k, v in my_ind_dict.items()})

    empty = GenomicPositionData.__new__(GenomicPositionData)
    empty.files_set = [path.name]
    empty.file_names = [path.stem]
    empty.config = Config(
        promoter_start=start,
        promoter_end=stop,
        chromosome=chrom,
        offset=start,
        viz_backend=viz,
    )
    empty.individual_data = ind_df
    empty.allelic_data = empty._generate_chisquared_test_df()
    empty.positions = empty.individual_data.columns.tolist()
    empty.means = df[["start", "mean"]].T.rename(columns=df["start"]).drop("start")
    empty.means.columns = empty.means.columns.astype(str)
    empty.modes = df[["start", "mode"]].T.rename(columns=df["start"]).drop("start")
    empty.modes.columns = empty.modes.columns.astype(str)
    empty.diffs = df[["start", "diff"]].T.rename(columns=df["start"]).drop("start")
    empty.diffs.columns = empty.diffs.columns.astype(str)
    return empty
