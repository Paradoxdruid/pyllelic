#!/usr/bin/env python3

"""Configuration options for pyllelic."""

from pathlib import Path

base_directory: Path = Path.cwd()
promoter_file: Path = base_directory / "promoter.txt"
results_directory: Path = base_directory / "results"
bam_directory: Path = base_directory / "bam_output"
analysis_directory: Path = base_directory / "test"
promoter_start: str = "1293000"
promoter_end: str = "1296000"
chromosome: str = "5"
