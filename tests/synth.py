# SPDX-License-Identifier: MIT
"""Deterministic synthetic reference and assemblies for the fast tests (tests/data/*.rr3 were encoded from these).
Only random.Random(seed).random() is used: its sequence is fixed across Python versions.

    python3 tests/synth.py <dir>     writes synth_ref.fa, synthA.fa, synthB.fa
"""
import os
import random
import sys

COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def _rng(seed):
    r = random.Random(seed)
    return lambda k: int(r.random() * k)


def _rand_seq(ri, n, alphabet=b"ACGT"):
    return bytes(alphabet[ri(len(alphabet))] for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


def _mutate(ri, s, snp=300, indel=2000):
    out, i, n = bytearray(), 0, len(s)
    while i < n:
        e = ri(snp * indel)
        if e < indel:                                     # substitution
            out.append(b"ACGT"[(b"ACGT".find(s[i:i + 1]) + 1 + ri(3)) % 4] if s[i:i + 1] in (b"A", b"C", b"G", b"T") else s[i])
            i += 1
        elif e < indel + snp // 2:                        # insertion
            out += _rand_seq(ri, 1 + ri(5))
        elif e < indel + snp:                             # deletion
            i += 1 + ri(5)
        else:
            out.append(s[i])
            i += 1
    return bytes(out)


def _lower(ri, s, runs, maxlen):
    b = bytearray(s)
    for _ in range(runs):
        p = ri(len(b))
        l = 1 + ri(maxlen)
        b[p:p + l] = bytes(b[p:p + l]).lower()
    return bytes(b)


def _fasta(recs, lw):
    o = bytearray()
    for (h, s), w in zip(recs, lw):
        o += b">" + h + b"\n"
        for i in range(0, len(s), w):
            o += s[i:i + w] + b"\n"
    return bytes(o)


def reference():
    ri = _rng(20261004)
    a = _rand_seq(ri, 400_000)
    b = _rand_seq(ri, 300_000)
    c = _rand_seq(ri, 20_000) + b"N" * 3000 + _rand_seq(ri, 27_000)
    return [(b"chrA synthetic", a), (b"chrB", b), (b"chrC", c)]


def assembly(seed):
    ri = _rng(seed)
    (_, a), (_, b), (_, c) = reference()
    recs = []
    recs.append((b"ctg1 copy of chrA", _lower(ri, _mutate(ri, a[10_000:210_000]), 40, 3000)))
    recs.append((b"ctg2\treverse complement of chrB", _mutate(ri, rc(b[5_000:150_000]))))
    recs.append((b"ctg3 inversion", _mutate(ri, a[250_000:300_000] + rc(a[300_000:320_000]) + a[320_000:380_000])))
    unit = _rand_seq(ri, 171)
    novel = _rand_seq(ri, 30_000) + _mutate(ri, unit * 60, snp=50) + b"N" * 5000 + _rand_seq(ri, 4000)
    recs.append((b"ctg4 novel and tandem", _lower(ri, novel, 25, 800)))
    for k, n in enumerate((5, 100, 1023, 1024, 1025, 4097, 16385)):
        recs.append((b"short%d" % k, _mutate(ri, c[1000 * k:1000 * k + n])[:n]))
    iupac = bytearray(_mutate(ri, b[200_000:230_000]))
    for _ in range(50):
        iupac[ri(len(iupac))] = b"RYKMSWBDHVn"[ri(11)]
    recs.append((b"ctg5 iupac", bytes(iupac)))
    recs.append((b"ctg6 mixed", _lower(ri, c[15_000:30_000], 10, 2000)))
    return recs


def write(d):
    os.makedirs(d, exist_ok=True)
    ref = reference()
    files = {"synth_ref.fa": _fasta(ref, [80] * len(ref))}
    for nm, seed in (("synthA", 1), ("synthB", 2)):
        recs = assembly(seed)
        files[nm + ".fa"] = _fasta(recs, [60 if i % 2 == 0 else 70 for i in range(len(recs))])
    for f, data in files.items():
        with open(os.path.join(d, f), "wb") as o:
            o.write(data)
    return {f: os.path.join(d, f) for f in files}


if __name__ == "__main__":
    write(sys.argv[1])
