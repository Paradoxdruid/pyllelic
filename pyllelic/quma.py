#!/usr/bin/env python3
"""Tools to quantify methylation in reduced representation
   bisulfite sequencing reads."""

# Adapted from QUMA CLI: http://quma.cdb.riken.jp/
# Licensed under GPLv3
# initial Perl conversion by: https://freelancer.com/u/Zubayerskd
# refactoring by Paradoxdruid

import re
from dataclasses import dataclass
from io import StringIO
from typing import Dict, List, Tuple, Optional
from Bio.Align import substitution_matrices
from Bio import Align

# import logging

# logging.basicConfig(filename="quma_test.log", level=logging.DEBUG)

MAX_LINE_LENGTH: int = 60

MAT: str = """    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U
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
MATRIX: substitution_matrices.Array = substitution_matrices.read(StringIO(MAT))

ALPHABET: str = "ACGTURYMWSKDHBVNacgturymwskdhbvn"

# ##################################################################


@dataclass
class Result:
    """Dataclass of quma aligment comparison results."""

    qAli: str = ""
    gAli: str = ""
    val: str = ""
    perc: float = 0.0
    pconv: float = 0.0
    gap: int = 0
    menum: int = 0
    unconv: int = 0
    conv: int = 0
    match: int = 0
    aliMis: int = 0
    aliLen: int = 0


@dataclass
class Fasta:
    """Dataclass to wrap fasta results."""

    com: str = ""
    pos: Optional[str] = None
    seq: str = ""


@dataclass
class Reference:
    """Dataclass of quma analysis intermediates.

    Includes fasta sequence, quma results, directon of read, genomic direction,
    and whether result meets exclusion criteria.
    """

    fasta: Fasta
    res: Result
    dir: int
    gdir: int
    exc: int


class Quma:
    """Quma methylation analysis parser for bisulfite conversion DNA sequencing."""

    def __init__(self, gfile_contents: str, qfile_contents: str) -> None:
        self._gfile_contents: str = gfile_contents
        """str: Genomic alignment sequence."""

        self._qfile_contents: str = qfile_contents
        """str: query alignment sequence."""

        self._gseq: str = self._parse_genome()
        self._qseq: List[Fasta] = self._parse_biseq()

        _gfilepF: str = self._fasta_make(self._gseq, "genomeF")

        self._raw_data: List[Reference] = self._process_fasta_output(
            self._qseq, "queryF", "queryR", _gfilepF
        )

        self.values: str = self._format_output(self._gseq, self._raw_data)
        """QUMA output values in tabular form."""

    def _parse_genome(self) -> str:
        """Parse genome file, removing white spaces and extra returns.

        Returns:
            str: parsed and curated string of genome sequence.
        """

        seq: str = self._gfile_contents

        seq = re.sub(r"^[\r\s]+", "", seq)
        seq = re.sub(r"[\r\s]+$", "", seq)
        seq = re.sub(r"(\r|\n|\r\n){2}", "\r|\n|\r\n", seq)

        return self._parse_seq(seq)

    @staticmethod
    def _scrub_whitespace(string: str) -> str:
        """Remove whitespace and repeated newlines.

        Args:
            string (str): input fasta string.

        Returns:
            str: scrubbed fasta string.
        """
        return (
            string.strip()
            .replace("\r\n\r\n", "\r\n")
            .replace("\n\n", "\n")
            .replace("\r\r", "\r")
            .replace("\r\n", "\n")
            .replace("\r", "\n")
        )

    def _parse_biseq(self) -> List[Fasta]:
        """Parse bisulfite sequencing fasta file.

        Returns:
            List[Fasta]: list of Fasta objects of sequence reads
        """

        multi: str = self._qfile_contents
        multi = self._scrub_whitespace(multi)
        biseq: List[Fasta] = []

        for line in multi.splitlines():

            if ">" in line:
                fa: Fasta = Fasta()
                fa.com = line
                fa.com = re.sub(r"^>", "", fa.com)
                fa.com = re.sub(r"\s*$", "", fa.com)
                biseq.append(fa)
            else:
                line = self._check_char_in_allowed(line, ALPHABET)
                if line == "":
                    continue
                fa.seq += line.upper()

        return biseq

    def _parse_seq(self, seq: str) -> str:
        """Extract sequence strings from the string of a text file.

        Args:
            seq (str): stringified text file.

        Returns:
            str: sequence string
        """

        seq = seq.replace("\r\n", "\n").replace("\r", "\n").upper()

        reg = re.compile(r"^\s*>.*?\n", re.MULTILINE)
        if re.findall(reg, seq):
            seq = re.sub(r"^\s*>(.*)?\n", "", seq)

        elif re.findall(
            r"^ORIGIN\s*\n((\s+(\d+(?:\s+\w+)+))+)\s*\n//", seq, re.MULTILINE
        ):  # pragma: no cover lgtm [py/redos]
            seq = re.findall(reg, seq, re.MULTILINE)[0][0]

        elif re.findall(
            r"^SQ\s+SEQUENCE.*\n((\s+(?:\w+\s+)+\d+)+)\n\/\/",
            seq,
            re.MULTILINE,
        ):  # pragma: no cover lgtm [py/redos]
            seq = re.findall(
                r"^SQ\s+SEQUENCE.*\n((\s+(?:\w+\s+)+\d+)+)\n\/\/",
                seq,
                re.MULTILINE,
            )[0][
                0
            ]  # lgtm [py/redos]

        elif re.findall(
            r"\.\.\s*\n((\s+(\d+(?:\s+\w+)+))+)\s*", seq, re.MULTILINE
        ):  # pragma: no cover lgtm [py/redos]
            seq = re.findall(
                r"\.\.\s*\n((\s+(\d+(?:\s+\w+)+))+)\s*", seq, re.MULTILINE
            )[
                0
            ]  # lgtm [py/redos]
        elif re.findall(r"^\s*>.+\s.+", seq, re.MULTILINE):  # pragma: no cover
            seq = re.findall(r"^\s*>(.+?)\s(?=.+)", seq, re.MULTILINE)[0][0]
            _ = seq

        return self._check_char_in_allowed(seq, ALPHABET)

    @staticmethod
    def _check_char_in_allowed(seq: str, pattern: str) -> str:
        """Return only charcters in string present in pattern.

        Args:
            seq (str): sequence string
            pattern (str): string of allowed characters

        Returns:
            str: string with unallowed characters removed.
        """

        return "".join([each for each in seq if each in pattern])

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
        qseq: List[Fasta],
        qfileF: str,
        qfileR: str,
        gfilepF: str,
    ) -> List[Reference]:
        """Process fasta alignment.

        Args:
            qseq (List[Fasta]): Fasta of query sequence
            qfileF (str): query sequence forward read
            qfileR (str): query sequence reverse complement
            gfilepF (str): genomic sequence forward read

        Returns:
            List[Reference]: list of quma References
        """

        UNCONV: int = 5
        PCONV: float = 95.0
        MIS: int = 10
        PERC: float = 90.0

        data: List[Reference] = []
        pos: int = 0
        for fa in qseq:
            pos += 1
            fa.pos = str(pos)

            qfilepF: str = self._fasta_make(fa.seq, qfileF)
            qfilepR: str = self._fasta_make(self._rev_comp(fa.seq), qfileR)

            fwd_result: Result = self._align_seq_and_generate_stats(qfilepF, gfilepF)
            rev_result: Result = self._align_seq_and_generate_stats(qfilepR, gfilepF)

            result: Result
            final_direction: int
            result, final_direction = self._find_best_dataset(fwd_result, rev_result)
            genome_direction: int = 1

            ref: Reference = Reference(
                fasta=fa, res=result, dir=final_direction, gdir=genome_direction, exc=0
            )
            # Flag for exclusion
            if (
                result.unconv > UNCONV
                or result.pconv > PCONV
                or result.aliMis > MIS
                or result.perc > PERC
            ):
                ref.exc = 1

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

        temp: List[str] = list(seq)
        temp.reverse()
        seq = "".join(temp)

        mappings: Dict[str, str] = {
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

        new: str = ""
        for each in seq:
            try:
                new += mappings[each]
            except KeyError:  # pragma: no cover
                new += each

        return new

    @staticmethod
    def _matching_substrings(alignment: Align.PairwiseAlignment) -> Tuple[str, str]:
        """Find pairwise alignment substrings.

        Args:
            alignment (Align.PairwiseAlignment): pairwise alignment

        Returns:
            Tuple[str, str]: query and genomic aligned substrings
        """

        matches: str = str(alignment).splitlines()[1]
        q_matches: str = str(alignment).splitlines()[0]
        g_matches: str = str(alignment).splitlines()[2]

        left_start_index: int = len(matches) - len(matches.lstrip())
        right_end_index: int = len(matches) - len(matches.rstrip())

        q_substring: str
        g_substring: str
        if right_end_index == 0:
            q_substring = q_matches[left_start_index:]
            g_substring = g_matches[left_start_index:]
        else:
            q_substring = q_matches[left_start_index:-right_end_index]
            g_substring = g_matches[left_start_index:-right_end_index]

        return (q_substring, g_substring)

    def _align_seq_and_generate_stats(self, gfile: str, qfile: str) -> Result:
        """Run pairwise sequence alignment.

        Args:
            gfile (str): genomic sequence file contents
            qfile (str): sequencing read(s) file contents

        Returns:
            Result: results dataclass
        """

        result: Result = Result()

        bio_gseq: str = gfile.splitlines()[1]
        bio_qseq: str = qfile.splitlines()[1]

        aligner: Align.PairwiseAligner = Align.PairwiseAligner(
            mode="local",
            substitution_matrix=MATRIX,
            open_gap_score=-10,
            extend_gap_score=-0.5,
        )
        bio_alignments: Align.PairwiseAlignments = list(
            aligner.align(bio_qseq, bio_gseq)
        )
        bio_alignment: Align.PairwiseAlignment = bio_alignments[0]

        query_ali: str
        genome_ali: str
        query_ali, genome_ali = self._matching_substrings(bio_alignment)

        fh_: str = f">genome\n{genome_ali}\n>que\n{query_ali}\n"
        # logging.debug(f"gseq={bio_gseq}\nqseq={bio_qseq}\nOut:{fh_}")

        fh: List[str] = fh_.split("\n")

        for i, line in enumerate(fh):
            if ">que" in line:
                result.qAli = fh[i + 1].strip()
            elif ">genome" in line:
                result.gAli = fh[i + 1].strip()

        result.qAli = result.qAli.replace(" ", "-")
        result.gAli = result.gAli.replace(" ", "-")

        results: Result = self._process_alignment_matches(result)

        return results

    def _process_alignment_matches(self, result: Result) -> Result:
        """Process alignment data to populate results dictionary.

        Args:
            ref (Dict[str, Any]): initial results dictionary

        Returns:
            Dict[str, Any]: results dictionary
        """

        gAli: str = result.gAli
        qAli: str = result.qAli[: len(gAli)]

        cpg: List[int] = [m.start() for m in re.finditer("CG", gAli)]
        result.aliLen = len(qAli)

        # Loop through sequence looking for CpG conversion
        j: int = 0
        for i in range(result.aliLen):

            g: str = gAli[i]
            q: str = qAli[i]

            if g != "-":
                j += 1

            if q == g or g == "C" and q == "T":
                result.match += 1

            if g == "-":
                result.gap += 1
                continue

            if q == "-":
                result.gap += 1
                if (j - 1) in cpg:
                    result.val += "-"
                continue

            if i == result.aliLen - 1:
                break

            # End if genome sequence isn't a C
            if g != "C":
                continue

            # If this is a CpG C on the genome, inspect query
            if (j - 1) in cpg:
                if q == "C":
                    result.conv += 1
                    result.menum += 1
                    result.val += "1"
                elif q == "T":
                    result.unconv += 1
                    result.val += "0"
                else:
                    result.val += q
                continue

        # kludge:
        if result.val == "":
            result.val = "-"

        results: Result = self._generate_summary_stats(result)
        # logging.debug(f"results={results}")
        return results

    def _generate_summary_stats(self, result: Result) -> Result:
        """Helper to generate summary statistics in results dictionary.

        Args:
            result (Result): quma results

        Returns:
            Result: quma results with sumary statistics
        """

        if result.conv + result.unconv != 0:
            result.pconv = self._percentage(result.conv, result.unconv, calc_type="sum")
        else:
            result.pconv = 0

        result.perc = self._percentage(result.match, result.aliLen, calc_type="total")
        result.aliMis = result.aliLen - result.match

        return result

    @staticmethod
    def _percentage(a: int, b: int, calc_type: str) -> float:
        """Helper to return percentages.

        Args:
            a (int): numerator
            b (int): denominator
            calc_type (str): 'sum' or 'total' calculation type

        Returns:
           float: percentage
        """

        if calc_type == "sum":
            return float(f"{(100 * a / (a + b)):3.1f}")
        if calc_type == "total":
            return float(f"{(100 * a / b):3.1f}")
        raise ValueError("No or incorrect calc_type provided")

    @staticmethod
    def _find_best_dataset(ffres: Result, frres: Result) -> Tuple[Result, int]:
        """Helper to find best data returned.

        Args:
            ffres (Result): quma result from forward alignment
            frres (Result): quma result from reverse alignment

        Returns:
            Tuple[Result, int]: best quma result and direction
        """

        # Find best dataset: FIXME

        fres: Result
        fdir: int
        if ffres.aliLen > frres.aliLen:
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
    def _format_output(gseq: str, data: List[Reference]) -> str:
        """Process program output into quma-formatted string.

        Args:
            gseq (str): genomic sequence
            data (List[Reference]): list of quma result References

        Returns:
            str: output tabular quma-format string
        """

        CONV_OUTPUT: int = 0
        output: str = ""

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

        for reference in data:

            output += (
                "\t".join(
                    [
                        str(reference.fasta.pos),
                        str(reference.fasta.com),
                        str(reference.fasta.seq),
                        str(reference.res.qAli),
                        str(reference.res.gAli),
                        str(reference.res.aliLen),
                        str(reference.res.aliMis),
                        str(reference.res.perc),
                        str(reference.res.gap),
                        str(reference.res.menum),
                        str(reference.res.unconv),
                        str(reference.res.conv),
                        str(reference.res.pconv),
                        str(reference.res.val),
                        str(reference.dir),
                        str(reference.gdir),
                    ]
                )
                + "\n"
            )
        return output
