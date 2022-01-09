#!/usr/bin/env python3
"""Configuration options for pyllelic."""

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass
class Config:
    base_directory: Path = Path("/")
    promoter_file: Path = base_directory / "promoter.txt"
    results_directory: Path = base_directory / "results"
    analysis_directory: Path = base_directory / "test"
    promoter_start: int = 1293200
    promoter_end: int = 1296000
    chromosome: str = "5"
    offset: int = 1298163
    viz_backend: str = "plotly"
    fname_pattern: str = re.compile(r"^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$")
