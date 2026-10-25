# SPDX-License-Identifier: MIT
"""Fetch and verify a panvram q16k cohort from the Zenodo record, for the window-loader measurement.

The record (10.5281/zenodo.23232317) ships the 558 archives as one tar cut into 500 MiB parts; the individual
`*.q16k.rr3` files are NOT addressable on their own (they return 404). So:

  n_asm == 558: every part, checked against PARTS.sha256, joined, the whole tar checked against its .sha256,
                extracted, every archive checked against SHA256SUMS -> 558 archives.
  n_asm == 5:   the first 100 MB of part00 by HTTP range, the tar headers walked by hand, the archives that are
                complete inside that prefix extracted, the first five kept and checked against SHA256SUMS.q16k.

Either way: an archive of 0 bytes, a missing archive or a wrong SHA-256 raises RuntimeError at once.

Note for the five-assembly mode: the tar is ordered `./`, `./SHA256SUMS`, then the archives by name, so the five are
the first five of the release-2 assemblies (r2_HG00097_hap1, ...), NOT the five year-1 haplotypes of the CPU rows of
logs/loader-2026-10-22. The cohort size is what this mode controls; the names differ, and the report must say so.
"""
import hashlib
import os
import shutil
import subprocess
import tarfile
from pathlib import Path

REC = "https://zenodo.org/records/23232317/files"
PREFIX_BYTES = 100 * 1024 * 1024
BLOCK = 512


def _md5(path, buf=8 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _sha256(path, buf=8 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _fetch(url, dest, byte_range=None, tries=5):
    dest = Path(dest)
    for attempt in range(1, tries + 1):
        cmd = ["curl", "-fsS", "--retry", "3", "--retry-delay", "5", "-o", str(dest)]
        if byte_range:
            cmd += ["-r", byte_range]
        r = subprocess.run(cmd + [url], capture_output=True, text=True)
        if r.returncode == 0 and dest.exists() and dest.stat().st_size > 0:
            return dest
        if attempt == tries:
            raise RuntimeError(f"download failed after {tries} attempts: {url}\n"
                               f"curl said: {r.stderr.strip() or r.returncode}; "
                               f"size on disk: {dest.stat().st_size if dest.exists() else 'no file'}")
    raise AssertionError


def _read_sums(path):
    sums = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        digest, name = line.split(None, 1)
        sums[os.path.basename(name.strip())] = digest
    return sums


def _verify(dirpath, sums, expect_n):
    dirpath = Path(dirpath)
    got = sorted(p for p in dirpath.glob("*.q16k.rr3"))
    if len(got) != expect_n:
        raise RuntimeError(f"expected {expect_n} archives in {dirpath}, found {len(got)}")
    for p in got:
        if p.stat().st_size == 0:
            raise RuntimeError(f"{p.name} is 0 bytes - the download did not deliver the file")
        want = sums.get(p.name)
        if want is None:
            raise RuntimeError(f"{p.name} is not listed in the checksum file")
        have = _sha256(p)
        if have != want:
            raise RuntimeError(f"{p.name}: SHA-256 {have} does not match the record's {want}")
    print(f"VERIFIED\t{len(got)} archives, every SHA-256 equals the record's, no empty file")
    for p in got[:3]:
        print(f"  {p.name}\t{p.stat().st_size} B\tsha256 {_sha256(p)}\tmd5 {_md5(p)}")
    if len(got) > 3:
        print(f"  ... {len(got) - 3} more")
    return got


def _walk_tar_prefix(buf_path, out_dir, want):
    """Extract the archives wholly contained in a truncated tar, in tar order, up to `want` of them."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    taken, members = [], []
    with open(buf_path, "rb") as f:
        size = Path(buf_path).stat().st_size
        off = 0
        while off + BLOCK <= size and len(taken) < want:
            f.seek(off)
            head = f.read(BLOCK)
            if head[:1] == b"\0":
                break
            name = head[0:100].rstrip(b"\0").decode("utf-8", "replace")
            raw = head[124:136].rstrip(b"\0 ").decode() or "0"
            msize = int(raw, 8)
            data_off = off + BLOCK
            padded = (msize + BLOCK - 1) // BLOCK * BLOCK
            members.append((name, msize))
            base = os.path.basename(name)
            if base.endswith(".q16k.rr3"):
                if data_off + msize > size:
                    break                                    # truncated member: stop, do not write a short file
                f.seek(data_off)
                (out_dir / base).write_bytes(f.read(msize))
                taken.append(base)
            elif base == "SHA256SUMS" and data_off + msize <= size:
                f.seek(data_off)
                (out_dir / "SHA256SUMS.in_tar").write_bytes(f.read(msize))
            off = data_off + padded
    print(f"PREFIX\t{len(members)} tar members seen, {len(taken)} archives complete: {', '.join(taken)}")
    if len(taken) < want:
        raise RuntimeError(f"only {len(taken)} complete archives in the first {PREFIX_BYTES} bytes, need {want}; "
                           f"raise PREFIX_BYTES")
    return taken


def prepare_cohort(n_asm, dest="cohort", work="dl"):
    """Put a verified q16k cohort of n_asm assemblies in `dest`. Returns the list of archive paths."""
    if n_asm not in (5, 558):
        raise ValueError("n_asm must be 5 or 558")
    dest, work = Path(dest), Path(work)
    dest.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    if n_asm == 558:
        parts = [f"cohort558_q16k.tar.part{i:02d}" for i in range(16)]
        _fetch(f"{REC}/PARTS.sha256", work / "PARTS.sha256")
        psums = _read_sums(work / "PARTS.sha256")
        for p in parts:
            q = work / p
            if not (q.exists() and q.stat().st_size and _sha256(q) == psums.get(p)):
                _fetch(f"{REC}/{p}", q)
                have = _sha256(q)
                if have != psums.get(p):
                    raise RuntimeError(f"{p}: SHA-256 {have} does not match PARTS.sha256 {psums.get(p)}")
            print(f"  part ok {p} {q.stat().st_size} B", flush=True)
        tar = work / "cohort558_q16k.tar"
        with open(tar, "wb") as out:
            for p in parts:
                with open(work / p, "rb") as f:
                    shutil.copyfileobj(f, out, 8 << 20)
        _fetch(f"{REC}/cohort558_q16k.tar.sha256", work / "tar.sha256")
        want_tar = open(work / "tar.sha256").read().split()[0]
        have_tar = _sha256(tar)
        if have_tar != want_tar:
            raise RuntimeError(f"joined tar SHA-256 {have_tar} does not match the record's {want_tar}")
        print(f"TAR\t{tar.stat().st_size} B, SHA-256 matches the record", flush=True)
        with tarfile.open(tar) as t:
            t.extractall(dest, filter="data")
        for p in parts:
            (work / p).unlink(missing_ok=True)
        tar.unlink(missing_ok=True)
        inner = dest / "SHA256SUMS"
        if not inner.exists():
            raise RuntimeError("the tar did not contain SHA256SUMS")
        return _verify(dest, _read_sums(inner), 558)

    prefix = work / "part00.prefix"
    _fetch(f"{REC}/cohort558_q16k.tar.part00", prefix, byte_range=f"0-{PREFIX_BYTES - 1}")
    print(f"RANGE\t{prefix.stat().st_size} B of part00", flush=True)
    _walk_tar_prefix(prefix, dest, want=5)
    _fetch(f"{REC}/SHA256SUMS.q16k", work / "SHA256SUMS.q16k")
    out = _verify(dest, _read_sums(work / "SHA256SUMS.q16k"), 5)
    prefix.unlink(missing_ok=True)
    return out
