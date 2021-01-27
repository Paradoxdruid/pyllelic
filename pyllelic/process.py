#!/usr/bin/env python3
"""Utilities to pre-process and prepare data for use in pyllelic.
"""

import gzip
from Bio import SeqIO
from pathlib import Path
from typing import List, Dict, Any


def process_fastq_to_list(filepath: Path) -> List[Any]:
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

    record_list: List[Any] = []
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


def make_records_to_dictionary(record_list: List[Any]) -> Dict[str, Any]:
    """Take in list of biopython SeqRecords and output a dictionary
       with keys of the record name.

    Args:
        record_list (List[SeqRecord]): biopython sequence records from a fastq file

    Returns:
        Dict[str, SeqRecord]: dict of biopython SeqRecords from a fastq file
    """
    return dict(zip([record.id for record in record_list], record_list))
