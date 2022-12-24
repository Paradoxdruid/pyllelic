#!/usr/bin/env python3
"""Encoded inputs for testing pyllelic."""

# flake8: noqa

import numpy as np
import pandas as pd

# Test data
SAMPLE_DICT_OF_DFS = {
    "TEST1": pd.DataFrame(
        {
            "1": ["111", "111", "111", "011", "011", "011"],
            "2": ["111", "111", "111", "011", "011", "011"],
            "3": ["111", "111", "111", "111", "111", "111"],
        }
    ),
    "TEST2": pd.DataFrame(
        {
            "1": ["111", "111", "111", "111", "111", "111"],
            "2": ["111", "111", "111", "111", "111", "111"],
            "3": ["111", "111", "111", "111", "111", "111"],
        }
    ),
    "TEST3": pd.DataFrame(
        {
            "1": ["FAIL", "FAIL", "FAIL", "FAIL", "FAIL", "FAIL"],
            "2": ["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"],
            "3": ["111", "111", "111", "011", "011", "NaN"],
        }
    ),
}

TEST_DF_SEQ_READS = pd.DataFrame.from_dict(
    {
        "1": "ATGCGTACGTA",
        "2": "ATGCTAGGCCC",
    },
    orient="index",
    columns=["sequence"],
).stack()

EXPECTED_RAW_QUMA = (
    "genome\t0\tATCGTAGTCGA\t1\t0\n1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\t"
    + "ATCGTAGTCGA\t11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n2\tquery2\t"
    + "ATCGATAGCATT\tATCG-TAGT\tATCGATAGC\t9\t1\t88.9\t1\t1\t0\t1\t"
    + "100.0\t1\t1\t1\n"
)

EXPECTED_INTERMEDIATE_MEANS = pd.DataFrame.from_dict(
    {
        "TEST1": [np.float64(5 / 6), np.float64(5 / 6), np.float64(1.0)],
        "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
        "TEST3": [np.nan, np.nan, np.float64(13 / 15)],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_INTERMEDIATE_MODES = pd.DataFrame.from_dict(
    {
        "TEST1": [np.float64(2 / 3), np.float64(2 / 3), np.float64(1.0)],
        "TEST2": [np.float64(1.0), np.float64(1.0), np.float64(1.0)],
        "TEST3": [np.nan, np.nan, np.float64(1.0)],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_MEANS = pd.DataFrame.from_dict(
    {
        "1293588": {"test.bam": 1.0},
        "1294262": {"test.bam": 1.0},
        "1294316": {"test.bam": 1.0},
        "1294343": {"test.bam": 1.0},
        "1294369": {"test.bam": 1.0},
        "1294419": {"test.bam": 0.6666666666666666},
        "1294872": {"test.bam": 1.0},
        "1294972": {"test.bam": 1.0},
        "1295089": {"test.bam": 1.0},
        "1295246": {"test.bam": 1.0},
        "1295320": {"test.bam": 1.0},
        "1295365": {"test.bam": 1.0},
        "1295743": {"test.bam": 1.0},
        "1295770": {"test.bam": 1.0},
        "1295876": {"test.bam": 1.0},
        "1295903": {"test.bam": 1.0},
        "1295937": {"test.bam": 1.0},
        "1295964": {"test.bam": 1.0},
    }
)

EXPECTED_INTERMEDIATE_DIFFS = pd.DataFrame.from_dict(
    {
        "TEST1": [np.float64(1 / 6), np.float64(1 / 6), np.float64(0.0)],
        "TEST2": [np.float64(0.0), np.float64(0.0), np.float64(0.0)],
        "TEST3": [np.nan, np.nan, np.float64(-2 / 15)],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA = pd.DataFrame.from_dict(
    {
        "TEST1": [
            [1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)],
            [1.0, 1.0, 1.0, (2 / 3), (2 / 3), (2 / 3)],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ],
        "TEST2": [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ],
        "TEST3": [np.nan, np.nan, [1.0, 1.0, 1.0, (2 / 3), (2 / 3)]],
    },
    orient="index",
    columns=["1", "2", "3"],
)

EXPECTED_QUMA_VALUES = "genome\t0\tATCGTAGTCGA\t1\t0\n1\tquery1\tATCGTAGTCGA\tATCGTAGTCGA\tATCGTAGTCGA\t11\t0\t100.0\t0\t2\t0\t2\t100.0\t11\t1\t1\n2\tquery2\tATCGATAGCATT\tATCGATAGC\tATCG-TAGT\t9\t2\t77.8\t1\t1\t0\t1\t100.0\t1\t1\t1\n"


INPUT_READ_FILE = [
    ">read0\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read1\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read2\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read3\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read4\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read5\nGGCGTAGGTAGGTTCGTACGAAGTCGTA",
]

EXPECTED_BAM_OUTPUT_VALUES = {
    "1293589": ">read0\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read1\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read2\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read3\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read4\nCGGCGTAGGTAGGTTCGTACGAAGTCGTA\n>read5\nTGGCGTAGGTAGGTTCGTACGAAGTCGTA",
    "1294263": ">read0\nCGGTTTAGGGGTAGCGTTACGTTTGGGTT",
    "1294317": ">read0\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read1\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read2\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read3\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read4\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read5\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read6\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read7\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read8\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read9\nAACACTACCCCCGCGCCTCCTCGCACCCG\n>read10\nAACACTACCCCCGCGCCTCCTCGCACCCG",
    "1294344": ">read0\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read1\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read2\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read3\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read4\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read5\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read6\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read7\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read8\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read9\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG\n>read10\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG",
    "1294370": ">read0\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read1\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read2\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read3\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read4\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read5\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read6\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read7\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read8\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read9\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read10\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read11\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read12\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read13\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read14\nCGGTTTTTTTGACGTTATGGTTTTAGGTT\n>read15\nCGGTTTTTTTGACGTTATGGTTTTAGGTT",
    "1294420": ">read0\nCGAAATCCACTAACGTATAACGAAAACCG\n>read1\nCGAAATCCACTAACGTATAACGAAAACCG\n>read2\nCGAAATCCACTAACGTATAACGAAAACCG\n>read3\nCGAAATCCACTAACGTATAACGAAAACCG\n>read4\nCGAAATCCACTAACGTATAACGAAAACCG\n>read5\nCGAAATCCACTAACGTATAACGAAAACCG\n>read6\nCGAAATCCACTAACGTATAACGAAAACCG\n>read7\nCGAAATCCACTAACGTATAACGAAAACCG\n>read8\nCGAAATCCACTAACGTATAACGAAAACCG\n>read9\nCGAAATCCACTAACGTATAACGAAAACCG\n>read10\nCGAAATCCACTAACGTATAACGAAAACCG\n>read11\nCGAAATCCACTAACGTATAACGAAAACCG\n>read12\nCGAAATCCACTAACGTATAACGAAAACCG\n>read13\nCGAAATCCACTAACGTATAACGAAAACCG\n>read14\nCGAAATCCACTAACGTATAACGAAAACCG\n>read15\nCGAAATCCACTAACGTATAACGAAAACCG\n>read16\nCGAAATCCACTAACGTATAACGAAAACCG",
    "1294873": ">read0\nCGGGGAGGTTTATTTGGCGGAAGGAGGGG\n>read1\nCGGGGAGGTTTATTTGGCGGAAGGAGGGG",
    "1294973": ">read0\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read1\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read2\nNGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read3\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read4\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read5\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read6\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read7\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read8\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read9\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read10\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT\n>read11\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT",
    "1295090": ">read0\nAAAAACGCGCGACATCGCAAAAATAACCG",
    "1295247": ">read0\nCGGGAGGGGTTGGGACGGGGCGGGGTTCG\n>read1\nTGGGAGGGGTCGGGACGGGGCGGGGTTCG\n>read2\nTGGGAGGGGTTGGGATGGGGTGGGGTTTG\n>read3\nCGGGAGGGGTTGGGACGGGGCGGGGTTCG\n>read4\nCGGGAGGGGTTGGGACGGGGCGGGGTTCG",
    "1295321": ">read0\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read1\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read2\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read3\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read4\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read5\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read6\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read7\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read8\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read9\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read10\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read11\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read12\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read13\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read14\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read15\nTGGGTTTTTAGTTTTTTCGTTATGTGGGA\n>read16\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read17\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read18\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read19\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read20\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read21\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read22\nTGGGTTTTTAGTTTTTTCGTTATGTGGGA\n>read23\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read24\nTGGGTTTTTAGTTTTTTTGTTATGTGGGA\n>read25\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read26\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read27\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read28\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read29\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read30\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read31\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read32\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read33\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read34\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read35\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read36\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read37\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read38\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read39\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA\n>read40\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA",
    "1295366": ">read0\nTCTATACCCGCGAATCCACTAAAAACCCG\n>read1\nTCTATACCCGCGAATCCACTAAAAACCCG\n>read2\nTCTATACCCGCGAATCCACTAAAAACCCG",
    "1295744": ">read0\nCTCCCTAAACGAACACGCGAAACCTCCCG",
    "1295771": ">read0\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read1\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read2\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read3\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read4\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read5\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read6\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read7\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read8\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read9\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read10\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read11\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read12\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read13\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read14\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read15\nCGGAGTGTTTTTTTGTAATATTTTTTTGT\n>read16\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read17\nCGGAGTGTTTTTTTGTAATATTTTTTCGC\n>read18\nCGGAGTGTTTTTTTGTAATATTTTTTCGC",
    "1295877": ">read0\nACCACCCCAAATCTATTAATCACCCACCG",
    "1295904": ">read0\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read1\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read2\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read3\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read4\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read5\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read6\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read7\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read8\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA\n>read9\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA",
    "1295938": ">read0\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read1\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read2\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read3\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read4\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read5\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read6\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read7\nAACCAAACGCTCCTACTAACCGCGCACCG\n>read8\nAACCAAACGCTCCTACTAACCGCGCACCG",
    "1295965": ">read0\nCGGGTGTTTTATATTAGTTATAACGGTTT\n>read1\nCGGGTGTTTTATATTAGTTATAACGGTTT\n>read2\nCGGGTGTTTTATATTAGTTATAACGGTTT\n>read3\nCGGGTGTTTTATATTAGTTATAACGGTTT",
}

EXPECTED_WRITE_DF_OUTPUT = {
    "1293589": ">read0\nTGGCGTAGGTAGGTTCGTACGAAGTCGTA",
    "1294263": ">read0\nCGGTTTAGGGGTAGCGTTACGTTTGGGTT",
    "1294317": ">read0\nAACACTACCCCCGCGCCTCCTCGCACCCG",
    "1294344": ">read0\nTGGGGTTGGTAGGTTTAGGGGGATTTCGG",
    "1294370": ">read0\nCGGTTTTTTTGACGTTATGGTTTTAGGTT",
    "1294420": ">read0\nCGAAATCCACTAACGTATAACGAAAACCG",
    "1294873": ">read0\nCGGGGAGGTTTATTTGGCGGAAGGAGGGG",
    "1294973": ">read0\nTGGGTTTTCGTGTTGTATTAGTTGTTAGT",
    "1295090": ">read0\nAAAAACGCGCGACATCGCAAAAATAACCG",
    "1295247": ">read0\nCGGGAGGGGTTGGGACGGGGCGGGGTTCG",
    "1295321": ">read0\nCGGGTTTTTAGTTTTTTCGTTACGTGGGA",
    "1295366": ">read0\nTCTATACCCGCGAATCCACTAAAAACCCG",
    "1295744": ">read0\nCTCCCTAAACGAACACGCGAAACCTCCCG",
    "1295771": ">read0\nCGGAGTGTTTTTTTGTAATATTTTTTCGC",
    "1295877": ">read0\nACCACCCCAAATCTATTAATCACCCACCG",
    "1295904": ">read0\nCGGGGCGGTTTCGTCGAGAAAGGGTGGGA",
    "1295938": ">read0\nAACCAAACGCTCCTACTAACCGCGCACCG",
    "1295965": ">read0\nCGGGTGTTTTATATTAGTTATAACGGTTT",
}

EXPECTED_BAM_OUTPUT_POSITIONS = [
    "1293589",
    "1294263",
    "1294317",
    "1294344",
    "1294370",
    "1294420",
    "1294873",
    "1294973",
    "1295090",
    "1295247",
    "1295321",
    "1295366",
    "1295744",
    "1295771",
    "1295877",
    "1295904",
    "1295938",
    "1295965",
]

EXPECTED_BAM_INDIVIDUAL_POSITIONS = [
    "1293589",
    "1293592",
    "1293604",
    "1293608",
    "1293614",
    "1294263",
    "1294277",
    "1294282",
    "1294328",
    "1294330",
    "1294338",
    "1294344",
    "1294370",
    "1294386",
    "1294420",
    "1294434",
    "1294873",
    "1294890",
    "1294981",
    "1294983",
    "1295095",
    "1295097",
    "1295099",
    "1295105",
    "1295247",
    "1295257",
    "1295262",
    "1295267",
    "1295274",
    "1295321",
    "1295338",
    "1295343",
    "1295374",
    "1295376",
    "1295393",
    "1295753",
    "1295759",
    "1295761",
    "1295771",
    "1295797",
    "1295904",
    "1295909",
    "1295915",
    "1295918",
    "1295945",
    "1295958",
    "1295960",
    "1295965",
    "1295969",
    "1295988",
]

EXPECTED_ALLELIC_DATA = pd.DataFrame(
    {
        "cellLine": {
            0: "test",
            1: "test",
            2: "test",
            3: "test",
            4: "test",
            5: "test",
            6: "test",
            7: "test",
            8: "test",
            9: "test",
            10: "test",
            11: "test",
            12: "test",
            13: "test",
            14: "test",
            15: "test",
            16: "test",
        },
        "position": {
            0: "1294328",
            1: "1294330",
            2: "1294338",
            3: "1294344",
            4: "1294370",
            5: "1294386",
            6: "1294420",
            7: "1294434",
            8: "1294981",
            9: "1294983",
            10: "1295321",
            11: "1295338",
            12: "1295343",
            13: "1295771",
            14: "1295797",
            15: "1295904",
            16: "1295965",
        },
        "ad_stat": {
            0: 0.0009111188771537126,
            1: 0.0009111188771537126,
            2: 0.0009111188771537126,
            3: 0.0009111188771537126,
            4: 0.0009111188771537126,
            5: 6.334248366623988e-05,
            6: 3.737981840170153e-05,
            7: 3.737981840170153e-05,
            8: 0.0005320055051392492,
            9: 0.0005320055051392492,
            10: 7.074463098970701e-10,
            11: 1.1236418084450895e-09,
            12: 4.60092429500253e-08,
            13: 7.744216431044088e-06,
            14: 2.2090496998585475e-05,
            15: 0.0009111188771537126,
            16: 0.0003114909767673841,
        },
        "p_crit": {
            0: 1,
            1: 1,
            2: 1,
            3: 1,
            4: 1,
            5: 1,
            6: 1,
            7: 1,
            8: 1,
            9: 1,
            10: 1,
            11: 1,
            12: 1,
            13: 1,
            14: 1,
            15: 1,
            16: 1,
        },
        "raw": {
            0: 0.0009111188771537126,
            1: 0.0009111188771537126,
            2: 0.0009111188771537126,
            3: 0.0009111188771537126,
            4: 0.0009111188771537126,
            5: 6.334248366623988e-05,
            6: 3.737981840170153e-05,
            7: 3.737981840170153e-05,
            8: 0.0005320055051392492,
            9: 0.0005320055051392492,
            10: 7.074463098970701e-10,
            11: 1.1236418084450895e-09,
            12: 4.60092429500253e-08,
            13: 7.744216431044088e-06,
            14: 2.2090496998585475e-05,
            15: 0.0009111188771537126,
            16: 0.0003114909767673841,
        },
    }
)

EXPECTED_BAM_OUTPUT_GENOME_VALUES = {
    "1293589": ">genome1293589\ncggcgcaggcaggcccgcacgaagccgtac",
    "1294263": ">genome1294263\ncggctcaggggcagcgccacgcctgggcct",
    "1294317": ">genome1294317\nggcactgcccccgcgcctcctcgcacccgg",
    "1294344": ">genome1294344\ncggggctggcaggcccagggggaccccggc",
    "1294370": ">genome1294370\ncggcctccctgacgctatggttccaggccc",
    "1294420": ">genome1294420\ncggggtccactagcgtgtggcgggggccgg",
    "1294873": ">genome1294873\ncggggaggcccacctggcggaaggaggggg",
    "1294973": ">genome1294973\ncgggtccccgcgctgcaccagccgccagcc",
    "1295090": ">genome1295090\ngggagcgcgcggcatcgcgggggtggccgg",
    "1295247": ">genome1295247\ncgggaggggtcgggacggggcggggtccgc",
    "1295321": ">genome1295321\ncgggtccccagtccctccgccacgtgggaa",
    "1295366": ">genome1295366\ntctgtgcccgcgaatccactgggagcccgg",
    "1295744": ">genome1295744\nctccctggacgggcacgcgggacctcccgg",
    "1295771": ">genome1295771\ncggagtgcctccctgcaacacttccccgcg",
    "1295877": ">genome1295877\naccaccccaaatctgttaatcacccaccgg",
    "1295904": ">genome1295904\ncggggcggtcccgtcgagaaagggtgggaa",
    "1295938": ">genome1295938\nagccaggcgctcctgctggccgcgcaccgg",
    "1295965": ">genome1295965\ncgggcgcctcacaccagccacaacggcctt",
}

EXPECTED_STACKED_BAM = pd.Series(
    {
        ("1293589", "sequence"): "TGGCGTAGGTAGGTTCGTACGAAGTCGTA",
        ("1294263", "sequence"): "CGGTTTAGGGGTAGCGTTACGTTTGGGTT",
        ("1294317", "sequence"): "AACACTACCCCCGCGCCTCCTCGCACCCG",
        ("1294344", "sequence"): "TGGGGTTGGTAGGTTTAGGGGGATTTCGG",
        ("1294370", "sequence"): "CGGTTTTTTTGACGTTATGGTTTTAGGTT",
        ("1294420", "sequence"): "CGAAATCCACTAACGTATAACGAAAACCG",
        ("1294873", "sequence"): "CGGGGAGGTTTATTTGGCGGAAGGAGGGG",
        ("1294973", "sequence"): "TGGGTTTTCGTGTTGTATTAGTTGTTAGT",
        ("1295090", "sequence"): "AAAAACGCGCGACATCGCAAAAATAACCG",
        ("1295247", "sequence"): "CGGGAGGGGTTGGGACGGGGCGGGGTTCG",
        ("1295321", "sequence"): "CGGGTTTTTAGTTTTTTCGTTACGTGGGA",
        ("1295366", "sequence"): "TCTATACCCGCGAATCCACTAAAAACCCG",
        ("1295744", "sequence"): "CTCCCTAAACGAACACGCGAAACCTCCCG",
        ("1295771", "sequence"): "CGGAGTGTTTTTTTGTAATATTTTTTCGC",
        ("1295877", "sequence"): "ACCACCCCAAATCTATTAATCACCCACCG",
        ("1295904", "sequence"): "CGGGGCGGTTTCGTCGAGAAAGGGTGGGA",
        ("1295938", "sequence"): "AACCAAACGCTCCTACTAACCGCGCACCG",
        ("1295965", "sequence"): "CGGGTGTTTTATATTAGTTATAACGGTTT",
    }
)

SAMPLE_BED = """#chr	start	end	strand	type	total_num	methy_num	percent_num
chr1	10471	10472	+	CG	1	1	1.0
chr1	10476	10477	+	CHH	1	0	0.0
chr1	10477	10478	+	CHH	1	0	0.0
chr1	10478	10479	+	CHH	1	0	0.0
chr1	10480	10481	+	CHG	1	0	0.0
chr1	10483	10484	+	CHG	1	0	0.0
chr1	10484	10485	+	CG	1	1	1.0
chr1	10487	10488	+	CHH	1	0	0.0
chr1	10488	10489	+	CHG	1	0	0.0
"""

EXPECTED_METHBANK_DATA_MEANS = pd.DataFrame.from_dict(
    {
        "10471": {"mean": 1.0},
        "10476": {"mean": 0.0},
        "10477": {"mean": 0.0},
        "10478": {"mean": 0.0},
        "10480": {"mean": 0.0},
        "10483": {"mean": 0.0},
        "10484": {"mean": 1.0},
        "10487": {"mean": 0.0},
        "10488": {"mean": 0.0},
    }
)

# fmt: off
EXPECTED_INDIVIDUAL_POSITIONS = pd.DataFrame.from_dict(
    {1293589: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1293592: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1293604: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1293608: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1293614: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294263: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294277: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294282: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294328: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294330: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294338: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294344: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294370: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294386: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294420: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294434: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294873: {0: '1', 1: '1', 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294890: {0: '1', 1: '1', 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294981: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1294983: {0: '0', 1: '0', 2: '0', 3: '0', 4: '0', 5: '0', 6: '0', 7: '0', 8: '0', 9: '0', 10: '0', 11: '0', 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295095: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295097: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295099: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295105: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295247: {0: '1', 1: '1', 2: '1', 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295257: {0: '0', 1: '0', 2: '0', 3: '0', 4: '1', 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295262: {0: '0', 1: '1', 2: '1', 3: '1', 4: '1', 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295267: {0: '0', 1: '1', 2: '1', 3: '1', 4: '1', 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295274: {0: '1', 1: '1', 2: '1', 3: '1', 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295321: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: '1', 18: '1', 19: '1', 20: '1', 21: '1', 22: '1', 23: '1', 24: '1', 25: '1', 26: '1', 27: '1', 28: '1', 29: '1', 30: '1', 31: '1', 32: '1', 33: '1', 34: '1', 35: '1', 36: '1', 37: '1', 38: None, 39: None, 40: None}, 1295338: {0: '0', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: '1', 18: '1', 19: '1', 20: '1', 21: '1', 22: '1', 23: '1', 24: '1', 25: '1', 26: '1', 27: '1', 28: '1', 29: '1', 30: '1', 31: '1', 32: '1', 33: '1', 34: '1', 35: '1', 36: '1', 37: '1', 38: '1', 39: '1', 40: '1'}, 1295343: {0: '0', 1: '0', 2: '0', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: '1', 18: '1', 19: '1', 20: '1', 21: '1', 22: '1', 23: '1', 24: '1', 25: '1', 26: '1', 27: '1', 28: '1', 29: '1', 30: '1', 31: '1', 32: '1', 33: '1', 34: '1', 35: '1', 36: '1', 37: '1', 38: '1', 39: '1', 40: '1'}, 1295374: {0: '1', 1: '1', 2: '1', 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295376: {0: '1', 1: '1', 2: '1', 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295393: {0: '1', 1: '1', 2: '1', 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295753: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295759: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295761: {0: '1', 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295771: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: '1', 18: '1', 19: '1', 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295797: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: '1', 14: '1', 15: '1', 16: '1', 17: '1', 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295904: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295909: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295915: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295918: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295945: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295958: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295960: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295965: {0: '1', 1: '1', 2: '1', 3: '1', 4: '1', 5: '1', 6: '1', 7: '1', 8: '1', 9: '1', 10: '1', 11: '1', 12: '1', 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295969: {0: '0', 1: '0', 2: '0', 3: '0', 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}, 1295988: {0: '1', 1: '1', 2: '1', 3: '1', 4: None, 5: None, 6: None, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: None, 28: None, 29: None, 30: None, 31: None, 32: None, 33: None, 34: None, 35: None, 36: None, 37: None, 38: None, 39: None, 40: None}}
)
# fmt: on

TEST_PROM_FILE = """tcccacccccaggaagagggggttctcgtccccacctctcattcccacccttgaaattgcgaagaggattataggtaacctgcaggcaccctcgccagagcgtctgtgcttccagacacttctccccattgccggcaacccggctccactgccgcgcccagcctcctctgttcactgctctggcctcggcgcctggaaaccgcgtgtccatcaaaacgtgaaggtgaacctcgtaagtttatgcaaactggacaggagggagagcagaggcagagatcaccgtgtccactcgacgtcctgagcgaaaagccacgtgtgcccacgtgacgatggagacaggaggaccagggctctgcctgcccccttttctgagcccctactgcattcagctctggggcctgggccctcgacggccaccacctcctcacctgggctcctgcgcagccaagcgcagtcccgcacgctcatcttccacgtcagctcctgcagcgagagcttggcatgcttccccagggagatgaacttcttggtgttcctgaggaagcggcgttcgttgtgcctggagccccagaggcctgggggcaccagccggcgcaggcaggcccgcacgaagccgtacacctgccaggggctgctgtgctggcggagcagctgcaccaggcgacgggggtctgtgtcctcctcctcgggggccgccacagagccctggggcttctcccgggcacagacaccggctgctggggtgaccgcagctcgcagcgggcagtgcgtcttgaggagcaccccgtaggggcactgcgcgtggttcccaagcagctccagaaacaggggccgcatttgccagtagcgctggggcaggcggggcaacctgcggggagtccctggcatccagggcctggaacccagaaagatggtctccacgagcctccgagcgccagtcaggctgggcctcagagagctgagtaggaaggagggccgcagctgctccttgtcgcctgaggagtagaggaagtgcttggtctcggcgtacaccgggggacaaggcgtgtcccagggacgtggtggccgcgatgtggatggggggcccgcgtggtgctggcggcccacggatgggtgggagtggcgcgtgccagagagcgcaccctccaaagaggtggcttcttcggcgggtctggcaggtgacaccacacagaaaccacggtcactcggtccacgcgtcctgcccgggtgggcccaggacccctgcccaacgggcgtccgctccggctcaggggcagcgccacgcctgggcctcttgggcaacggcagacttcggctggcactgcccccgcgcctcctcgcacccggggctggcaggcccagggggaccccggcctccctgacgctatggttccaggcccgttcgcatcccagacgccttcggggtccactagcgtgtggcgggggccgggcctgagtggcagcgccgagctggtacagcggcggcccgcacacctggtaggcgcagctgggagccaccagcacaaagagcgcgcagcgtgccagcaggtgaaccagcacgtcgtcgcccacgcggcgcagcagcagcccccacgccccgctcccccgcagtgcgtcggtcaccgtgttgggcaggtagctgcgcacgctggtggtgaaggcctcgggggggcccccgcgggccccgtccagcagcgcgaagccgaaggccagcacgttcttcgcgccgcgctcgcacagcctctgcagcactcgggccaccagctccttcaggcaggacacctgcgggggaagcgccctgagtcgcctgcgctgctctccgcatgtcgctggttccccccggccgccctcaaccccagccggacgccgaccccggggaggcccacctggcggaaggagggggcggcggggggcggccgtgcgtcccagggcacgcacaccaggcactgggccaccagcgcgcggaaagccgccgggtccccgcgctgcaccagccgccagccctggggccccaggcgccgcacgaacgtggccagcggcagcacctcgcggtagtggctgcgcagcagggagcgcacggctcggcagcggggagcgcgcggcatcgcgggggtggccggggccagggcttcccacgtgcgcagcaggacgcagcgctgcctgaaactcgcgccgcgaggagagggcggggccgcggaaaggaaggggaggggctgggagggcccggagggggctgggccggggacccgggaggggtcgggacggggcggggtccgcgcggaggaggcggagctggaaggtgaaggggcaggacgggtgcccgggtccccagtccctccgccacgtgggaagcgcggtcctgggcgtctgtgcccgcgaatccactgggagcccggcctggccccgacagcgcagctgctccgggcggacccgggggtctgggccgcgcttccccgcccgcgcgccgctcgcgctcccagggtgcagggacgccagcgagggccccagcggagagaggtcgaatcggcctaggctgtggggtaacccgagggaggggccatgatgtggaggccctgggaacaggtgcgtgcggcgaccctttggccgctggcctgatccggagacccagggctgcctccaggtccggacgcggggcgtcgggctccgggcaccacgaatgccggacgtgaaggggaggacggaggcgcgtagacgcggctggggacgaacccgaggacgcattgctccctggacgggcacgcgggacctcccggagtgcctccctgcaacacttccccgcgacttgggctccttgacacaggcccgtcatttctctttgcaggttctcaggcggcgaggggtccccaccatgagcaaaccaccccaaatctgttaatcacccaccggggcggtcccgtcgagaaagggtgggaaatggagccaggcgctcctgctggccgcgcaccgggcgcctcacaccagccacaacggccttgaccct"""

SAMPLE_BAI = b"QkFJAVUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAUgIAAAEAAAAAAMwFAAAAAAAAfmsAAAAASpIAAAIAAAAAAMwFAAAAAAAAfmsAAAAAJgMAAAAAAAAAAAAAAAAAAFAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAAAAzAUAAAAAAADMBQAAAAAAAMwFAAAAAL1zdTEAAAAAAAAAAAAAAAA="

SAMPLE_BAM = b"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"
