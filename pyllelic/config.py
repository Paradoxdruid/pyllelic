#!/usr/bin/env python3
"""Configuration options for pyllelic."""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class Config:
    base_directory: Path = Path("/")
    promoter_file: Path = base_directory / "promoter.txt"
    results_directory: Path = base_directory / "results"
    bam_directory: Path = base_directory / "bam_output"
    analysis_directory: Path = base_directory / "test"
    promoter_start: str = "1293200"
    promoter_end: str = "1296000"
    chromosome: str = "5"
    offset: int = 1298163
