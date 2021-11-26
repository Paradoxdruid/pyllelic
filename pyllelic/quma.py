#!/usr/bin/env python3
"""Tools to quantify methylation in reduced representation
   bisulfite sequencing reads."""

# Adapted from QUMA CLI: http://quma.cdb.riken.jp/
# Licensed under GPLv3
# perl conversion by: https://freelancer.com/u/Zubayerskd

import re
from io import StringIO
from typing import Any, Dict, List, Tuple, Optional

from Bio import pairwise2
from Bio.Align import substitution_matrices

# import logging

# logging.basicConfig(filename="quma_test.log", level=logging.DEBUG)

MAX_LINE_LENGTH = 60
MAT = """    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U
    A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4
    T  -4   5  -4   5  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5
    G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4
    C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4
    S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4
    W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1   1
    R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4
    Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1   1
    K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1   1
    M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4
    B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -1
    V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4
    H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -1
    D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -1
    N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
    U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5
    """
MATRIX = substitution_matrices.read(StringIO(MAT))

ALPHABET = "ACGTURYMWSKDHBVNacgturymwskdhbvn"

# ##################################################################


class Quma:
    """Quma methylation analysis parser for bisulfite conversion DNA sequencing."""

    def __init__(self, gfile_contents: str, qfile_contents: str) -> None:
        self._gfile_contents = gfile_contents
        """str: Genomic alignment sequence."""

        self._qfile_contents = qfile_contents
        """str: query alignment sequence."""

        self._gseq = self._parse_genome()
        self._qseq = self._parse_biseq()

        cpgf = self._find_cpg(self._gseq)

        gfilepF = self._fasta_make(self._gseq, "genomeF")

        self._raw_data: List[Dict[str, Any]] = self._process_fasta_output(
            self._qseq, "queryF", "queryR", gfilepF, cpgf
        )

        self.values = self._format_output(self._gseq, self._raw_data)
        """QUMA output values in tabular form."""

    def _parse_genome(self) -> str:
        """Parse genome file, removing white spaces and extra returns.

        Returns:
            str: parsed and curated string of genome sequence.
        """

        seq = self._gfile_contents

        seq = re.sub(r"^[\r\s]+", "", seq)
        seq = re.sub(r"[\r\s]+$", "", seq)
        seq = re.sub(r"(\r|\n|\r\n){2}", "\r|\n|\r\n", seq)

        return self._parse_seq(seq)

    def _parse_biseq(self) -> List[Dict[str, str]]:
        """Parse bisulfite sequencing fasta file.

        Returns:
            List[Dict[str, str]]: list of dictionaries of sequence reads
        """

        multi = self._qfile_contents
        multi = re.sub(r"^[\r\s]+", "", multi)
        multi = re.sub(r"[\r\s]+$", "", multi)
        multi = re.sub(r"(\r\n){2}", "\r\n", multi)
        multi = re.sub(r"(\n){2}", "\n", multi)
        multi = re.sub(r"(\r){2}", "\r", multi)

        multi = re.sub(r"\r\n", "\n", multi)
        multi = re.sub(r"\r", "\n", multi)
        biseq: List[Dict[str, str]] = []
        fa: Dict[str, str] = {}

        multi = re.findall(r"(.*)$", multi, re.MULTILINE)

        for line in multi:

            if ">" in line:
                if fa and not fa["seq"]:  # pragma: no cover
                    biseq.pop()
                fa = {"com": line}
                fa["com"] = re.sub(r"^>", "", fa["com"])
                fa["com"] = re.sub(r"\s*$", "", fa["com"])
                biseq.append(fa)
            else:
                line = self._check_char_in_allowed(line, ALPHABET)
                if line == "":
                    continue
                # if not fa:  # pragma: no cover
                #     return None  # does this ever happen?
                try:
                    fa["seq"] += line.upper()
                except KeyError:
                    fa["seq"] = line.upper()

        if fa:
            if not fa.get("seq"):  # pragma: no cover
                biseq.pop()

        return biseq

    def _parse_seq(self, seq: str) -> str:
        """Extract sequence strings from the string of a text file.

        Args:
            seq (str): stringified text file.

        Returns:
            str: sequence string
        """
        _ = ""
        seq = re.sub(r"\r\n", "\n", seq)
        seq = re.sub(r"\r", "\n", seq)
        seq = seq.upper()

        reg = re.compile(r"^\s*>.*?\n", re.MULTILINE)
        if re.findall(reg, seq):
            seq = re.sub(r"^\s*>(.*)?\n", "", seq)

        elif re.findall(
            r"^ORIGIN\s*\n((\s+(\d+(\s+\w+)+))+)\s*\n//", seq, re.MULTILINE
        ):  # pragma: no cover
            seq = re.findall(reg, seq, re.MULTILINE)[0]

        elif re.findall(
            r"^SQ\s+SEQUENCE.*\n((\s+(\w+\s+)+\d+)+)\n\/\/", seq, re.MULTILINE
        ):  # pragma: no cover
            seq = re.findall(
                r"^SQ\s+SEQUENCE.*\n((\s+(\w+\s+)+\d+)+)\n\/\/", seq, re.MULTILINE
            )[0]

        elif re.findall(
            r"\.\.\s*\n((\s+(\d+(\s+\w+)+))+)\s*", seq, re.MULTILINE
        ):  # pragma: no cover
            seq = re.findall(r"\.\.\s*\n((\s+(\d+(\s+\w+)+))+)\s*", seq, re.MULTILINE)[
                0
            ]
        elif re.findall(r"^\s*>.+\s.+", seq, re.MULTILINE):  # pragma: no cover
            seq = re.findall(r"^\s*>(.+?)\s(?=.+)", seq, re.MULTILINE)[0]
            _ = seq

        return self._check_char_in_allowed(seq, ALPHABET)

    @staticmethod
    def _check_char_in_allowed(seq: str, pattern: str) -> str:
        """Return only charcters in string present in pattern.

        Args:
            seq (str): sequence string
            patt (str): string of allowed characters

        Returns:
            str: string with unallowed characters removed.
        """
        new = ""
        for each in seq:
            if each in pattern:
                new += each
        return new

    @staticmethod
    def _find_cpg(gseq: str) -> Tuple[List[str], Dict[str, Any], Dict[str, Any]]:
        """Find CpG sites in genomic string."""
        positions = []
        cpgf = {}
        cpgr = {}
        length = len(gseq)
        pos = 0
        while True:
            try:
                pos = gseq.index("CG", pos)
            except ValueError:
                break

            cpgf[str(pos)] = 1
            positions.append(str(pos))
            cpgr[str(length - pos - 2)] = 1
            pos += 1

        return cpgf

    @staticmethod
    def _fasta_make(seq: str, seq_name: str, line: int = None) -> str:
        """Write a sequence string to a fasta-formatted text file contents.

        Args:
            seq (str): sequence string
            seq_name (str): sequence name
            line (int, optional): Max line length to process. Defaults to None.

        Returns:
            str: fasta-formatted text file contents.
        """
        line = line or MAX_LINE_LENGTH

        seq = re.sub(r"[0-9]| |\t|\n|\r|\f", "", seq)

        reg = r"(.{1," + str(line) + "})"
        seq = re.sub(reg, r"\1\n", seq)

        return f">{seq_name}\n{seq}"

    def _process_fasta_output(
        self,
        qseq: List[Dict[str, str]],
        qfileF: str,
        qfileR: str,
        gfilepF: str,
        cpgf: Dict[str, int],
    ) -> List[Dict[str, Any]]:
        """Process fasta alignment."""

        UNCONV = 5
        PCONV = 95.0
        MIS = 10
        PERC = 90.0

        data: List[Dict[str, Any]] = []
        pos = 0
        for fa in qseq:
            pos += 1
            fa["pos"] = str(pos)

            qfilepF = self._fasta_make(fa["seq"], qfileF)
            qfilepR = self._fasta_make(self._rev_comp(fa["seq"]), qfileR)

            fwd_result = self._align_seq_and_generate_stats(qfilepF, gfilepF, cpgf)
            rev_result = self._align_seq_and_generate_stats(qfilepR, gfilepF, cpgf)

            result, final_direction = self._find_best_dataset(fwd_result, rev_result)
            genome_direction = 1

            ref: Dict[str, Any] = {
                "fa": fa,
                "res": result,
                "dir": final_direction,
                "gdir": genome_direction,
                "exc": 0,
            }
            # Flag for exclusion
            if (
                result["unconv"] > UNCONV
                or result["pconv"] > PCONV
                or result["aliMis"] > MIS
                or result["perc"] > PERC
            ):
                ref["exc"] = 1

            data.append(ref)
        return data

    @staticmethod
    def _rev_comp(seq: str) -> str:
        """Return reverse complement of sequence.

        Args:
            seq (str): sequence

        Returns:
            str: reverse complement of sequence
        """
        temp = list(seq)
        temp.reverse()
        seq = "".join(temp)

        mappings = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "U": "A",
            "R": "Y",
            "Y": "R",
            "M": "K",
            "W": "W",
            "S": "S",
            "K": "M",
            "D": "H",
            "H": "D",
            "B": "V",
            "V": "B",
            "N": "N",
            "a": "t",
            "c": "g",
            "g": "c",
            "t": "a",
            "u": "a",
            "r": "y",
            "y": "r",
            "m": "k",
            "w": "w",
            "s": "s",
            "k": "m",
            "d": "h",
            "h": "d",
            "b": "v",
            "v": "b",
            "n": "n",
        }

        new = ""
        for each in seq:
            try:
                new += mappings[each]
            except KeyError:  # pragma: no cover
                new += each

        return new

    def _align_seq_and_generate_stats(
        self, gfile: str, qfile: str, cpg: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Run pairwise sequence alignment.

        Args:
            gfile (str): genomic sequence file contents
            qfile (str): sequencing read(s) file contents
            cpg (Dict[str, Any]): dictionary of CpG locations

        Returns:
            Dict[str, Any]: results dictionary
        """

        ref: Dict[str, Any] = {
            "qAli": "",
            "gAli": "",
            "gap": 0,
            "menum": 0,
            "unconv": 0,
            "conv": 0,
            "pconv": "",
            "match": 0,
            "val": "",
            "perc": "",
            "aliMis": 0,
        }

        bio_gseq = gfile.splitlines()[1]
        bio_qseq = qfile.splitlines()[1]
        # flipping to correct output comparison
        bio_alignments = pairwise2.align.localds(bio_qseq, bio_gseq, MATRIX, -10, -0.5)
        genome_ali = bio_alignments[0].seqB[
            bio_alignments[0].start : bio_alignments[0].end
        ]
        query_ali = bio_alignments[0].seqA[
            bio_alignments[0].start : bio_alignments[0].end
        ]
        fh_ = f">genome\n{genome_ali}\n>que\n{query_ali}\n"
        # logging.debug(f"gseq={bio_gseq},aligned={genome_ali}")
        fh: List[str] = fh_.split("\n")

        for i, line in enumerate(fh):
            if ">que" in line:
                ref["qAli"] = fh[i + 1].strip()
            elif ">genome" in line:
                ref["gAli"] = fh[i + 1].strip()

        ref["qAli"] = ref["qAli"].replace(" ", "-")
        ref["gAli"] = ref["gAli"].replace(" ", "-")

        if len(ref["qAli"]) != len(ref["gAli"]):  # pragma: no cover
            print("qAli len != gAli len")
            # sys.exit()

        results = self._process_alignment_matches(ref, cpg)
        return results

    def _process_alignment_matches(
        self, ref: Dict[str, Any], cpg: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Process alignment data to populate results dictionary.

        Args:
            ref (Dict[str, Any]): initial results dictionary
            cpg (Dict[str, Any]): dictionary of CpG locations

        Returns:
            Dict[str, Any]: results dictionary
        """
        gAli = ref["gAli"]
        qAli = ref["qAli"][: len(gAli)]

        cpg = self._find_cpg(gAli)
        ref["aliLen"] = len(qAli)

        # Loop through sequence looking for CpG conversion
        j = 0
        for i in range(ref["aliLen"]):

            g = gAli[i]
            q = qAli[i]

            if g != "-":
                j += 1

            if q == g or g == "C" and q == "T":
                ref["match"] += 1

            if g == "-":
                ref["gap"] += 1
                continue

            if q == "-":
                ref["gap"] += 1
                if cpg.get(str(j - 1)):
                    ref["val"] += "-"
                continue

            if i == ref["aliLen"] - 1:
                break

            # End if genome sequence isn't a C
            if g != "C":
                continue

            # If this is a CpG C on the genome, inspect query
            if cpg.get(str(j - 1)):
                if q == "C":
                    ref["conv"] += 1
                    ref["menum"] += 1
                    ref["val"] += "1"
                elif q == "T":
                    ref["unconv"] += 1
                    ref["val"] += "0"
                else:
                    ref["val"] += q
                continue

        results = self._generate_summary_stats(ref)

        # kludge:
        if results["val"] == "":
            results["val"] = "-"
        # logging.debug(f"results={results}")
        return results

    def _generate_summary_stats(self, ref: Dict[str, Any]) -> Dict[str, Any]:
        """Helper to generate summary statistics in results dictionary."""

        if ref["conv"] + ref["unconv"] != 0:
            ref["pconv"] = self._percentage(ref["conv"], ref["unconv"], calc_type="sum")
        else:
            ref["pconv"] = 0

        ref["perc"] = self._percentage(ref["match"], ref["aliLen"], calc_type="total")
        ref["perc"] = float(ref["perc"])
        ref["pconv"] = float(ref["pconv"])
        ref["aliMis"] = ref["aliLen"] - ref["match"]

        return ref

    @staticmethod
    def _percentage(a: int, b: int, calc_type: str) -> Optional[str]:
        """Helper to return percentages."""
        if calc_type == "sum":
            return f"{(100 * a / (a + b)):3.1f}"
        if calc_type == "total":
            return f"{(100 * a / b):3.1f}"
        return None  # pragma: no cover

    @staticmethod
    def _find_best_dataset(
        ffres: Dict[str, Any], frres: Dict[str, Any]
    ) -> Tuple[Dict[str, Any], int]:
        """Helper to find best data returned."""
        # Find best dataset: FIXME

        if ffres["aliLen"] > frres["aliLen"]:
            fres = ffres
            fdir = 1
        else:
            fres = frres
            fdir = -1
        # print(f"Forward:\n{ffres}\n\nRevese:\n{frres}\n")
        # if ffres["aliMis"] > frres["aliMis"]:
        #     print("Cond 1")
        #     fres = frres
        #     fdir = -1
        # elif ffres["aliMis"] < frres["aliMis"]:
        #     print("Cond 2")
        #     fres = ffres
        #     fdir = 1
        # elif ffres["perc"] > frres["perc"]:
        #     print("Cond 3")
        #     fres = ffres
        #     fdir = 1
        # elif ffres["perc"] < frres["perc"]:
        #     print("Cond 4")
        #     fres = frres
        #     fdir = -1
        # elif ffres["unconv"] > frres["unconv"]:
        #     print("Cond 5")
        #     fres = frres
        #     fdir = -1
        # elif ffres["unconv"] < frres["unconv"]:
        #     print("Cond 6")
        #     fres = ffres
        #     fdir = 1
        # elif ffres["pconv"] < frres["pconv"]:
        #     print("Cond 7")
        #     fres = frres
        #     fdir = -1
        # elif ffres["pconv"] > frres["pconv"]:
        #     print("Cond 8")
        #     fres = ffres
        #     fdir = 1
        # else:
        #     print("Cond 9")
        #     fres = ffres
        #     fdir = 1

        # print(f"fres: {fres},\nfdir: {fdir}")
        return fres, fdir

    @staticmethod
    def _format_output(gseq, data) -> str:
        """Process program output into quma-formatted string."""
        CONV_OUTPUT = 0
        output = ""

        output += (
            "\t".join(
                [
                    "genome",
                    str(CONV_OUTPUT),
                    gseq,
                    "1",
                    ",".join(["0"]),
                ]
            )
            + "\n"
        )

        for ref in data:

            output += (
                "\t".join(
                    [
                        str(ref["fa"]["pos"]),
                        ref["fa"]["com"],
                        ref["fa"]["seq"],
                        ref["res"]["qAli"],
                        ref["res"]["gAli"],
                        str(ref["res"]["aliLen"]),
                        str(ref["res"]["aliMis"]),
                        str(ref["res"]["perc"]),
                        str(ref["res"]["gap"]),
                        str(ref["res"]["menum"]),
                        str(ref["res"]["unconv"]),
                        str(ref["res"]["conv"]),
                        str(ref["res"]["pconv"]),
                        str(ref["res"]["val"]),
                        str(ref["dir"]),
                        str(ref["gdir"]),
                    ]
                )
                + "\n"
            )
        return output
