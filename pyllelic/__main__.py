#!/usr/bin/env python3
"""pyllelic: module level interface to run pyllelic from the command line.

Example usage:

    python -m pyllelic -o my_data -f fh_cellline_tissue.fastq.gz -g hg19chr5 \
                        -chr chr5 -s 1293000 -e 1296000 --viz plotly

This command would save pyllelic results in files with the prefix `my_data`, analyzing
the specified fastq file using the specified reference genome, in the genomic region
indicated.
"""

from . import pyllelic
from . import process
from .config import Config
from pathlib import Path
import argparse
from typing import List


def _parsing() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: parsed arguments.
    """

    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="run pyllelic on bisulfite sequencing fastq files",
        prog="python -m pyllelic",
    )

    parser.add_argument("-o", "--output_fname", type=str, required=True)
    parser.add_argument("-f", "--fastq", type=str, required=True)
    parser.add_argument("-g", "--genome", type=str, required=True)
    parser.add_argument("-chr", "--chrom", type=str, required=True)
    parser.add_argument("-s", "--start", type=int, required=True)
    parser.add_argument("-e", "--end", type=int, required=True)
    parser.add_argument("--viz", type=str, default="plotly")
    parser.add_argument(
        "--fname_pattern", type=str, default="^[a-zA-Z]+_([a-zA-Z0-9]+)_.+bam$"
    )
    parser.add_argument("--testdir", type=str, default="test")

    args: argparse.Namespace = parser.parse_args()
    return args


def _process_files(args: argparse.Namespace) -> None:
    """Process fastq and genome files using Bismark.

    Args:
        args (argparse.Namespace): parsed args
    """

    process.retrieve_promoter_seq(
        "genome.txt", chrom=args.chrom, start=args.start, end=args.end
    )
    genome: Path = Path.cwd() / args.genome
    fastq: Path = Path.cwd() / args.fastq
    process.prepare_genome(genome)
    process.bismark(genome, fastq)

    bamfile: Path = Path.cwd() / (Path(args.fastq).stem + ".bam")
    process.pysam_sort(bamfile)
    process.pysam_index(Path(bamfile.parent) / (bamfile.stem + "_sorted.bam"))


def _call_pyllelic(args: argparse.Namespace) -> pyllelic.GenomicPositionData:
    """Run pyllelic data analysis.

    Args:
        args (argparse.Namespace): parsed args

    Returns:
        GenomicPositionData: pyllelic data object
    """

    fname_pattern: str = fr"{args.fname_pattern}"

    config: Config = pyllelic.configure(
        base_path=str(Path.cwd()),
        prom_file="genome.txt",
        prom_start=args.start,
        prom_end=args.end,
        chrom=args.chrom,
        offset=args.start,
        viz_backend=args.viz,
        fname_pattern=fname_pattern,
        test_dir=args.testdir,
    )

    files_set: List[str] = pyllelic.make_list_of_bam_files(config)

    data: pyllelic.GenomicPositionData = pyllelic.pyllelic(
        config=config, files_set=files_set
    )
    return data


def run_pyllelic() -> None:
    """Run all processing and analysis steps of pyllelic."""

    args: argparse.Namespace = _parsing()
    print("Preparing genome and processing fastq file... this will take a while.")
    _process_files(args)
    print("Running pyllelic...")
    data: pyllelic.GenomicPositionData = _call_pyllelic(args)
    print("Data processed, saving output files...")
    data.save_pickle(args.output_fname + ".pickle")
    data.save(args.output_fname + ".xlsx")
    print("Pyllelic run complete")


# Run the whole process
if __name__ == "__main__":  # pragma: no cover
    run_pyllelic()
