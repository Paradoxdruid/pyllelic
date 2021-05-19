# Adapted from QUMA CLI: http://quma.cdb.riken.jp/
# Licensed under GPLv3
# perl conversion by: https://freelancer.com/u/Zubayerskd

import time
from datetime import datetime
import re
import os
from io import StringIO
from typing import List, Tuple, Dict, Any
from Bio import pairwise2
from Bio.Align import substitution_matrices

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


def check_char_in_allowed(seq: str, pattern: str) -> str:
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


def make_time(given_time: float = None) -> str:
    """Convert floating time into string formatted timestamp.

    Args:
        given_time (float, optional): floating point time. Defaults to None.

    Returns:
        str: String-formatted timestamp.
    """
    curr_time: float = given_time or time.time()

    dt: datetime = datetime.fromtimestamp(curr_time)
    return dt.strftime("%Y%m%d%H%M%S")


def curate_seq(seq: str) -> str:
    """Curate a sequence to only have allowed characters.

    Args:
        seq (str): sequence to check.

    Returns:
        str: sequence without dis-allowed characters.
    """
    return check_char_in_allowed(seq, "ACGTURYMWSKDHBVNacgturymwskdhbvn")


def parse_seq(seq: str) -> str:
    """Extract sequence strings from the string of a text file.

    Args:
        seq (str): stringified text file.

    Returns:
        str: sequence string
    """
    _ = ""  # was com
    seq = re.sub(r"\r\n", "\n", seq)
    seq = re.sub(r"\r", "\n", seq)
    seq = seq.upper()

    reg = re.compile(r"^\s*>.*?\n", re.MULTILINE)
    if re.findall(reg, seq):
        seq = re.sub(r"^\s*>(.*)?\n", "", seq)

    elif re.findall(r"^ORIGIN\s*\n((\s+(\d+(\s+\w+)+))+)\s*\n//", seq, re.MULTILINE):
        seq = re.findall(reg, seq, re.MULTILINE)[0]

    elif re.findall(r"^SQ\s+SEQUENCE.*\n((\s+(\w+\s+)+\d+)+)\n\/\/", seq, re.MULTILINE):
        seq = re.findall(
            r"^SQ\s+SEQUENCE.*\n((\s+(\w+\s+)+\d+)+)\n\/\/", seq, re.MULTILINE
        )[0]

    elif re.findall(r"\.\.\s*\n((\s+(\d+(\s+\w+)+))+)\s*", seq, re.MULTILINE):
        seq = re.findall(r"\.\.\s*\n((\s+(\d+(\s+\w+)+))+)\s*", seq, re.MULTILINE)[0]
    elif re.findall(r"^\s*>.+\s.+", seq, re.MULTILINE):
        seq = re.findall(r"^\s*>(.+?)\s(?=.+)", seq, re.MULTILINE)[0]
        _ = seq  # was com

    return curate_seq(seq)


def parse_genome(file: str) -> str:
    """Parse genome file, removing white spaces and extra returns.

    Args:
        file (str): text file path for fasta file.

    Returns:
        str: parsed and curated string of genome sequence.
    """
    with open(file, "r") as f:
        seq = f.read()

    seq = re.sub(r"^[\r\s]+", "", seq)
    seq = re.sub(r"[\r\s]+$", "", seq)
    seq = re.sub(r"(\r|\n|\r\n){2}", "\r|\n|\r\n", seq)

    return parse_seq(seq)


def multi_fasta_parse(multi: Any) -> List[Dict[str, str]]:
    """Find all bisufite sequencing reads in a fasta file,
       and return as a dictionary.

    Args:
        multi (Any): multi-line sequence string

    Returns:
        List[Dict[str, str]]: list of dictionaries of sequence reads
    """
    multi = re.sub(r"\r\n", "\n", multi)
    multi = re.sub(r"\r", "\n", multi)
    biseq: List[Dict[str, str]] = []
    fa: Dict[str, str] = {}

    multi = re.findall(r"(.*)$", multi, re.MULTILINE)

    for line in multi:

        if ">" in line:
            if fa and not fa["seq"]:
                biseq.pop()
            fa = {"com": line}
            fa["com"] = re.sub(r"^>", "", fa["com"])
            fa["com"] = re.sub(r"\s*$", "", fa["com"])
            biseq.append(fa)
        else:
            line = curate_seq(line)
            if line == "":
                continue
            if not fa:
                return None  # does this ever happen?
            try:
                fa["seq"] += line.upper()
            except KeyError:
                fa["seq"] = line.upper()

    if fa:
        try:
            _ = fa["seq"]  # was temp
        except KeyError:
            biseq.pop()
    return biseq


def parse_biseq(file: str) -> List[Dict[str, str]]:
    """Parse bisulfite sequencing fasta file.

    Args:
        file (str): file path.

    Returns:
        List[Dict[str, str]]: list of dictionaries of sequence reads
    """
    with open(file, "r") as f:
        multi = f.read()

    multi = re.sub(r"^[\r\s]+", "", multi)
    multi = re.sub(r"[\r\s]+$", "", multi)
    multi = re.sub(r"(\r\n){2}", "\r\n", multi)
    multi = re.sub(r"(\n){2}", "\n", multi)
    multi = re.sub(r"(\r){2}", "\r", multi)

    return multi_fasta_parse(multi)


def parse_multi(file: str) -> Tuple[None, List[Dict[str, str]]]:
    """Parse fast sequencing file with multiple reads.

    Args:
        file (str): file path

    Returns:
        Tuple[None, List[Dict[str, str]]]: None and list of dicts of sequence reads
    """
    with open(file, "r") as f:
        multi = f.read()

    multi = re.sub(r"^[\r\s]+", "", multi)
    multi = re.sub(r"[\r\s]+$", "", multi)
    multi = re.sub(r"(\r\n){2}", "\r\n", multi)
    multi = re.sub(r"(\n){2}", "\n", multi)
    multi = re.sub(r"(\r){2}", "\r", multi)

    biseq = multi_fasta_parse(multi)
    return None, biseq


def fasta_make(seq: str, seq_name: str, line: int = None) -> str:
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


def fasta_print(
    seq: str, seq_name: str, path: str, line: int = None, add: bool = None
) -> None:
    """Write a fast-formatted sequence file.

    Args:
        seq (str): sequence
        seq_name (str): sequence name
        path (str): path to write to
        line (int, optional): max length to process. Defaults to None.
        add (bool, optional): append to file flag. Defaults to None.
    """
    if add:
        with open(path, "a") as f:
            f.write(fasta_make(seq, seq_name, line))

    else:
        with open(path, "w") as f:
            f.write(fasta_make(seq, seq_name, line))


def fasta_output(seq: str, seq_name: str, line: int = None) -> str:
    """Return a fasta-formatted sequence file as a string

    Args:
        seq (str): sequence
        seq_name (str): sequence name
        line (int, optional): max length to process. Defaults to None.

    Returns:
        str: string of fasta-formatted file.
    """
    return fasta_make(seq, seq_name, line)


def rev_comp(seq: str) -> str:
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
        except KeyError:
            new += each

    return new


def align_seq_and_generate_stats(
    gfile: str, qfile: str, cpg: Dict[str, Any]
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
    bio_alignments = pairwise2.align.globalds(bio_qseq, bio_gseq, MATRIX, -10, -0.5)
    fh_ = f">genome\n{bio_alignments[0][0]}\n>que\n{bio_alignments[0][1]}\n"
    fh: List[str] = fh_.split("\n")

    for i, line in enumerate(fh):
        if ">que" in line:
            ref["qAli"] = fh[i + 1].strip()
        elif ">genome" in line:
            ref["gAli"] = fh[i + 1].strip()

    ref["qAli"] = ref["qAli"].replace(" ", "-")
    ref["gAli"] = ref["gAli"].replace(" ", "-")

    if len(ref["qAli"]) != len(ref["gAli"]):
        print("qAli len != gAli len")
        # sys.exit()

    results = process_alignment_matches(ref, cpg)
    return results


def process_alignment_matches(
    ref: Dict[str, Any], cpg: Dict[str, Any]
) -> Dict[str, Any]:
    """Process alignment data to populate results dictionary.

    Args:
        ref (Dict[str, Any]): initial results dictionary
        cpg (Dict[str, Any]): dictionary of CpG locations

    Returns:
        Dict[str, Any]: results dictionary
    """
    # ref = ref
    # cpg = cpg
    gAli = qAli = ""
    # fl = 0

    for i, val in enumerate(ref["gAli"]):
        # removed gap stripping, not needed for alignment function
        # if val == "-" and fl == 0:
        #     continue
        # fl = 1
        gAli += val
        qAli += ref["qAli"][i]

    # gAli = gAli.rstrip("-")
    qAli = qAli[: len(gAli)]

    ref["aliLen"] = len(qAli)

    j = 0
    for i in range(ref["aliLen"]):

        g = gAli[i]
        q = qAli[i]

        # if g != "-":
        #     j += 1
        j += 1
        if q == g or (g == "C" and q == "T"):
            ref["match"] += 1

        if g == "-":
            ref["gap"] += 1
            continue

        if q == "-":
            ref["gap"] += 1
            try:
                _ = cpg[str(j - 1)]
                ref["val"] += "-"
            except KeyError:
                pass
            continue

        if i == ref["aliLen"] - 1:
            break

        if g != "C":
            continue

        try:
            _ = cpg[str(j - 1)]
            if q == "C":
                ref["conv"] += 1
            elif q == "T":
                ref["unconv"] += 1
        except KeyError:
            continue

        if q == "C":
            ref["menum"] += 1
            ref["val"] += "1"
        elif q == "T":
            ref["val"] += "0"
        else:
            ref["val"] += q

    if ref["conv"] + ref["unconv"] != 0:
        ref["pconv"] = "{:3.1f}".format(
            100 * ref["conv"] / (ref["conv"] + ref["unconv"])
        )
    else:
        ref["pconv"] = 0

    ref["perc"] = "{:3.1f}".format(100 * ref["match"] / ref["aliLen"])
    ref["perc"] = float(ref["perc"])
    ref["pconv"] = float(ref["pconv"])

    ref["aliMis"] = ref["aliLen"] - ref["match"]
    return ref


def process_fasta_output(
    qseq: List[Dict[str, str]],
    qfileF: str,
    qfileR: str,
    gfilepF: str,
    gfilepR: str,
    cpgf: Dict[str, int],
    cpgr: Dict[str, int],
) -> List[Dict[str, Any]]:
    """Process fasta alignment."""

    UNCONVL = 5
    PCONVL = 95.0
    MISL = 10
    PERCL = 90.0

    # conv = 0
    unc = UNCONVL
    pcon = PCONVL
    mis = MISL
    perc = PERCL

    data: List[Dict[str, Any]] = []
    pos = 0
    for fa in qseq:
        pos += 1
        fa["pos"] = str(pos)  # was int

        qfilepF = fasta_output(fa["seq"], qfileF)
        qfilepR = fasta_output(rev_comp(fa["seq"]), qfileR)

        # if conv != 1:  # it always is
        ffres = align_seq_and_generate_stats(qfilepF, gfilepF, cpgf)
        frres = align_seq_and_generate_stats(qfilepR, gfilepF, cpgf)

        # Find best dataset:
        if ffres["aliMis"] > frres["aliMis"]:
            fres = frres
            fdir = -1
        elif ffres["aliMis"] < frres["aliMis"]:
            fres = ffres
            fdir = 1
        elif ffres["perc"] > frres["perc"]:
            fres = ffres
            fdir = 1
        elif ffres["perc"] < frres["perc"]:
            fres = frres
            fdir = -1
        elif ffres["unconv"] > frres["unconv"]:
            fres = frres
            fdir = -1
        elif ffres["unconv"] < frres["unconv"]:
            fres = ffres
            fdir = 1
        elif ffres["pconv"] < frres["pconv"]:
            fres = frres
            fdir = -1
        elif ffres["pconv"] > frres["pconv"]:
            fres = ffres
            fdir = 1
        else:
            fres = ffres
            fdir = 1

        res = fres
        dir = fdir
        gdir = 1

        ref: Dict[str, Any] = {"fa": fa, "res": res, "dir": dir, "gdir": gdir, "exc": 0}
        if res["unconv"] > unc:
            ref["exc"] = 1
        if res["pconv"] > pcon:
            ref["exc"] = 1
        if res["aliMis"] > mis:
            ref["exc"] = 1
        if res["perc"] > perc:
            ref["exc"] = 1

        data.append(ref)
    return data


def format_output(gseq, positions, data) -> str:
    """Process program output into quma-formatted string."""
    CONV_OUTPUT = 0
    output = ""

    output += (
        "\t".join(
            ["genome", str(CONV_OUTPUT), gseq, str(len(positions)), ",".join(positions)]
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


def find_cpg(gseq: str) -> Tuple[List[str], Dict[str, Any], Dict[str, Any]]:
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

    return (positions, cpgf, cpgr)


def quma_main(gfile: str, qfile: str) -> str:
    """Run quma for quantification of methylation of bisulfite sequencing reads."""

    t = make_time()
    uid = "{}{:06d}".format(t, os.getpid())

    gseq = parse_genome(gfile)
    qseq = parse_biseq(qfile)

    positions, cpgf, cpgr = find_cpg(gseq)

    qfileF = f"que{uid}F"
    qfileR = f"que{uid}R"
    gfileF = f"genome{uid}F"
    gfileR = f"genome{uid}R"

    gfilepF = fasta_output(gseq, gfileF)
    gfilepR = fasta_output(gseq, gfileR)

    data: List[Dict[str, Any]] = process_fasta_output(
        qseq, qfileF, qfileR, gfilepF, gfilepR, cpgf, cpgr
    )

    positions = list(map(str, positions))

    output = format_output(gseq, positions, data)

    return output
