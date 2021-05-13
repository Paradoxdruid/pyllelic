# Adapted from QUMA CLI: http://quma.cdb.riken.jp/
# Licensed under GPLv3
# Converted from perl and developed by: https://freelancer.com/u/Zubayerskd

import time
from datetime import datetime
import re
import os
import subprocess
import sys


LINE = 60


def y_cd(seq, patt):
    new = ""
    for each in seq:
        if each in patt:
            new += each
    return new


def usage():
    print(
        """
usage: quma.py [options] - or input_file or -g genome_file -q query_file

Input data
    1) - (STDIN)
        Multi-FASTA format of genomic DNA sequence and bisulfite sequences
        First sequence must be genomic DNA seuqnce
    2) input_file
        Multi-FASTA format of genomic DNA sequence and bisulfite sequences
        First sequence must be genomic DNA seuqnce
    3) -g genome_file -q query_file
        genome_file : FASTA format genomic DNA sequence
        query_file: Multi-FASTA format bisulfite sequences

Option
    -f: output format (0|1|2|3) default 0
        0: tab separated data
           first line : 'genome', condition of convert direction (see below),
                        genomic sequence, number of CpG, CpG position (first base = 0)
           other lines: No., sequence name, sequence, alignment data of this sequnece,
                        alignment data of genome sequence, alignment length,
                        number of alignmnet mismatch, percent identity of alignment,
                        number of alignment gap, number of methylated CpG,
                        number of bisulfite unconverted CpH (CpH: CpA, CpC, CpT),
                        number of bisulfite converted CpH,
                        percent of bisulfite convertion,
                        CpG methylation pattern (0: unmethylated, 1: methylated),
                        convert direction (1: forward, -1: reverse)
        1: human readable alignment data
        2: tab separated multiple alignment data
        3: tab separated summarized data
     -d: condition of convert direction of genomic sequence (0|1|2) default 0
        0: C -> T conversion
           PCR primer pair was designded for forward strand of the genomic sequence
        1: G -> A conversion
           PCR primer pair was designded for reverse strand of the genomic sequence
        2: both
           Search both direction of conversion and adopt more appropriate strand
     -u: upper limit of unconverted CpHs (integer, default 5)
         (CpH: CpA, CpC, CpT)
     -c: lower limit of percent converted CpHs (float, default 95.0)
     -m: upper limit of alignment mismatches (integer, default 10)
     -p: lower limit of percent identity (float, default 90.0)
     -u, -c, -m -p options are only affected output format 3 or 4"""
    )
    sys.exit()


def makeTime(given_time=None):
    curr_time = given_time or time.time()

    dt = datetime.fromtimestamp(curr_time)
    return dt.strftime("%Y%m%d%H%M%S")


def curateSeq(seq):
    return y_cd(seq, "ACGTURYMWSKDHBVNacgturymwskdhbvn")


def parseSeq(seq):
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

    return curateSeq(seq)


def parseGenome(file):
    with open(file, "r") as f:
        seq = f.read()

    seq = re.sub(r"^[\r\s]+", "", seq)
    seq = re.sub(r"[\r\s]+$", "", seq)
    seq = re.sub(r"(\r|\n|\r\n){2}", "\r|\n|\r\n", seq)

    return parseSeq(seq)


def multiFastaParse(multi):
    multi = re.sub(r"\r\n", "\n", multi)
    multi = re.sub(r"\r", "\n", multi)
    biseq = []
    fa = {}

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
            line = curateSeq(line)
            if line == "":
                continue
            if not fa:
                return
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


def parseBiseq(file):
    with open(file, "r") as f:
        multi = f.read()

    multi = re.sub(r"^[\r\s]+", "", multi)
    multi = re.sub(r"[\r\s]+$", "", multi)
    multi = re.sub(r"(\r\n){2}", "\r\n", multi)
    multi = re.sub(r"(\n){2}", "\n", multi)
    multi = re.sub(r"(\r){2}", "\r", multi)

    return multiFastaParse(multi)


def parseMulti(file):
    with open(file, "r") as f:
        multi = f.read()

    multi = re.sub(r"^[\r\s]+", "", multi)
    multi = re.sub(r"[\r\s]+$", "", multi)
    multi = re.sub(r"(\r\n){2}", "\r\n", multi)
    multi = re.sub(r"(\n){2}", "\n", multi)
    multi = re.sub(r"(\r){2}", "\r", multi)

    biseq = multiFastaParse(multi)
    return None, biseq


def fastaMake(seq, com, line=None):
    line = line or 60

    seq = re.sub(r"[0-9]| |\t|\n|\r|\f", "", seq)

    reg = r"(.{1," + str(line) + "})"
    seq = re.sub(reg, r"\1\n", seq)

    return f">{com}\n{seq}"


def fastaPrint(seq, com, path, line=None, add=None):
    if add:
        with open(path, "a") as f:
            f.write(fastaMake(seq, com, line))

    else:
        with open(path, "w") as f:
            f.write(fastaMake(seq, com, line))


def revComp(seq):
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


def execNeedle(gfile, qfile, cpg, needl, NDLOPT):

    ref = {
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

    # create the command for the neddle
    com = needl + f" {gfile} {qfile} " + NDLOPT

    fh = subprocess.run(com, shell=True, capture_output=True, text=True).stdout
    fh = fh.split("\n")

    for i, line in enumerate(fh):
        if ">que" in line:
            ref["qAli"] = fh[i + 1].strip()
        elif ">genome" in line:
            ref["gAli"] = fh[i + 1].strip()

    ref["qAli"] = ref["qAli"].replace(" ", "-")
    ref["gAli"] = ref["gAli"].replace(" ", "-")
    temp1 = ref["qAli"]
    temp2 = ref["gAli"]

    if len(temp1) != len(temp2):
        print("qAli len != gAli len")
        sys.exit()

    gAli = qAli = ""
    fl = 0

    for i, val in enumerate(ref["gAli"]):
        if val == "-" and fl == 0:
            continue
        fl = 1
        gAli += val
        qAli += ref["qAli"][i]

    gAli = gAli.rstrip("-")
    qAli = qAli[: len(gAli)]

    ref["aliLen"] = len(qAli)

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
            try:
                _ = cpg[str(j - 1)]  # was temp
                ref["val"] += "-"
            except KeyError:
                pass
            continue

        if i == ref["aliLen"] - 1:
            break

        if g != "C":
            continue

        try:
            _ = cpg[str(j - 1)]  # was temp
        except KeyError:
            if q == "C":
                ref["unconv"] += 1
            elif q == "T":
                ref["conv"] += 1
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


def quma_main(g_file, q_file, needle="needle", tempdir="/tmp/"):
    UNCONVL = 5
    PCONVL = 95.0
    MISL = 10
    PERCL = 90.0
    TEMPDIR = tempdir
    # NEEDLE = needle
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
    CPGMAT = TEMPDIR + "EDNA_CPG"
    NDLOPT = (
        " -stdout -auto -gapopen 10.0 -gapextend 0.5 "
        + "-aformat markx3 -datafile "
        + CPGMAT
    )
    t = makeTime()
    uid = "{}{:06d}".format(t, os.getpid())
    data = []
    positions = []

    file = None
    conv = None
    unc = None
    pcon = None
    mis = None
    perc = None
    gfile = g_file
    qfile = q_file

    conv = conv or 0
    unc = unc or UNCONVL
    pcon = pcon or PCONVL
    mis = mis or MISL
    perc = perc or PERCL

    if file:
        gseq, qseq = parseMulti(file)

    else:
        gseq = parseGenome(gfile)
        qseq = parseBiseq(qfile)

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
        positions.append(pos)
        cpgr[str(length - pos - 2)] = 1
        pos += 1

    with open(CPGMAT, "w") as f:
        f.write(MAT)

    qfileF = f"que{uid}F"
    qfileR = f"que{uid}R"
    gfileF = f"genome{uid}F"
    gfileR = f"genome{uid}R"
    qfilepF = TEMPDIR + qfileF
    qfilepR = TEMPDIR + qfileR
    gfilepF = TEMPDIR + gfileF
    gfilepR = TEMPDIR + gfileR

    fastaPrint(gseq, gfileF, gfilepF)
    fastaPrint(revComp(gseq), gfileR, gfilepR)

    pos = 0
    for fa in qseq:
        pos += 1
        fa["pos"] = pos

        fastaPrint(fa["seq"], qfileF, qfilepF)
        fastaPrint(revComp(fa["seq"]), qfileR, qfilepR)

        if conv != 1:
            ffres = execNeedle(qfilepF, gfilepF, cpgf, needle, NDLOPT)
            frres = execNeedle(qfilepR, gfilepF, cpgf, needle, NDLOPT)

            for _ in range(0, 1):
                if ffres["aliMis"] > frres["aliMis"]:
                    fres = frres
                    fdir = -1
                    break

                if ffres["aliMis"] < frres["aliMis"]:
                    fres = ffres
                    fdir = 1
                    break

                if ffres["perc"] > frres["perc"]:
                    fres = ffres
                    fdir = 1
                    break

                if ffres["perc"] < frres["perc"]:
                    fres = frres
                    fdir = -1
                    break

                if ffres["unconv"] > frres["unconv"]:
                    fres = frres
                    fdir = -1
                    break

                if ffres["unconv"] < frres["unconv"]:
                    fres = ffres
                    fdir = 1
                    break

                if ffres["pconv"] < frres["pconv"]:
                    fres = frres
                    fdir = -1
                    break

                if ffres["pconv"] > frres["pconv"]:
                    fres = ffres
                    fdir = 1
                    break

                fres = ffres
                fdir = 1

        if conv != 0:
            rfres = execNeedle(qfilepF, gfilepR, cpgr, needle, NDLOPT)
            rrres = execNeedle(qfilepR, gfilepR, cpgr, needle, NDLOPT)

            for t in range(0, 1):
                if rfres["aliMis"] > rrres["aliMis"]:
                    rres = rrres
                    rdir = -1
                    break

                if rfres["aliMis"] < rrres["aliMis"]:
                    rres = rfres
                    rdir = 1
                    break

                if rfres["perc"] > rrres["perc"]:
                    rres = rfres
                    rdir = 1
                    break

                if rfres["perc"] < rrres["perc"]:
                    rres = rrres
                    rdir = -1
                    break

                if rfres["unconv"] > rrres["unconv"]:
                    rres = rrres
                    rdir = -1
                    break

                if rfres["unconv"] < rrres["unconv"]:
                    rres = rfres
                    rdir = 1
                    break

                if rfres["pconv"] < rrres["pconv"]:
                    rres = rrres
                    rdir = -1
                    break

                if rfres["pconv"] > rrres["pconv"]:
                    rres = rfres
                    rdir = 1
                    break

            rres = rfres
            rdir = 1

        if conv == 0:
            res = fres
            dir = fdir
            gdir = 1

        elif conv == 1:
            res = rres
            dir = rdir
            gdir = -1

        else:
            if fres["aliMis"] > rres["aliMis"]:
                res = rres
                dir = rdir
                gdir = -1
                break

            if fres["aliMis"] < rres["aliMis"]:
                res = fres
                dir = fdir
                gdir = 1
                break

            if fres["perc"] > rres["perc"]:
                res = fres
                dir = fdir
                gdir = 1
                break

            if fres["perc"] < rres["perc"]:
                res = rres
                dir = rdir
                gdir = -1
                break

            if fres["unconv"] > rres["unconv"]:
                res = rres
                dir = rdir
                gdir = -1
                break

            if fres["unconv"] < rres["unconv"]:
                res = fres
                dir = fdir
                gdir = 1
                break

            if fres["pconv"] < rres["pconv"]:
                res = rres
                dir = rdir
                gdir = -1
                break

            if fres["pconv"] > rres["pconv"]:
                res = fres
                dir = fdir
                gdir = 1
                break

            res = fres
            dir = fdir
            gdir = 1

        if gdir != 1:
            res["qAli"] = revComp(res["qAli"])
            res["gAli"] = revComp(res["gAli"])
            temp = list(revComp(res["val"]))
            temp.reverse()
            res["val"] = "".join(temp)

        ref = {"fa": fa, "res": res, "dir": dir, "gdir": gdir, "exc": 0}
        if res["unconv"] > unc:
            ref["exc"] = 1
        if res["pconv"] > pcon:
            ref["exc"] = 1
        if res["aliMis"] > mis:
            ref["exc"] = 1
        if res["perc"] > perc:
            ref["exc"] = 1

        data.append(ref)

    os.unlink(gfilepF)
    os.unlink(gfilepR)
    os.unlink(qfilepF)
    os.unlink(qfilepR)
    os.unlink(CPGMAT)

    positions = list(map(str, positions))

    output = ""

    output += (
        "\t".join(["genome", str(conv), gseq, str(len(positions)), ",".join(positions)])
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
