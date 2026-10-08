# SPDX-License-Identifier: MIT
"""part_filter.py - `split --filter` program: receives one part on stdin ($FILE = part name), writes it to a temporary
file, records sha256 and md5 (parts table), uploads it to the Zenodo draft through the bucket PUT (curl, token via a
0600 config file; up to 5 attempts; a part already on the server with the same md5 is not re-sent), compares the md5
of the response with the local md5, deletes the temporary file. One part on disk at a time. Never publishes.
Env: ZENODO_STATE (deposition.json), PARTS_TABLE (tsv: name, bytes, sha256, md5, server md5, seconds), PART_DIR (tmp dir)."""
import hashlib, json, os, subprocess, sys, tempfile, time
from pathlib import Path
import urllib.request

name = os.environ["FILE"]; state = json.load(open(os.environ["ZENODO_STATE"])); table = Path(os.environ["PARTS_TABLE"])
TOKEN = Path("~/.zenodo_token").expanduser().read_text().strip(); d = Path(os.environ.get("PART_DIR", "/tmp"))
tmp = d / (name + ".tmp"); h256, hmd5, size = hashlib.sha256(), hashlib.md5(), 0
with open(tmp, "wb") as o:
    for b in iter(lambda: sys.stdin.buffer.read(1 << 24), b""): o.write(b); h256.update(b); hmd5.update(b); size += len(b)
sha, md5 = h256.hexdigest(), hmd5.hexdigest()


def server_files():
    """file list of the draft; 5 attempts (the listing itself answers 504 now and then)"""
    req = urllib.request.Request(f"https://zenodo.org/api/deposit/depositions/{state['id']}/files"); req.add_header("Authorization", "Bearer " + TOKEN)
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=120) as r: return {f["filename"]: f for f in json.load(r)}
        except Exception as e:
            print(f"{name}: file list attempt {attempt}: {type(e).__name__} {getattr(e, 'code', '')}", file=sys.stderr, flush=True); time.sleep(30 * attempt)
    return {}


t0 = time.time(); verdict = None
try:
    have = server_files().get(name)
    if have and have.get("checksum", "").replace("md5:", "") == md5:
        verdict = "already on server, md5 =="
    else:
        for attempt in range(1, 6):
            with tempfile.NamedTemporaryFile("w", delete=False, dir=str(d), prefix=".curlcfg", suffix=".tmp") as cfg:
                os.chmod(cfg.name, 0o600); cfg.write(f'header = "Authorization: Bearer {TOKEN}"\n')
            try:
                r = subprocess.run(["curl", "-sS", "-K", cfg.name, "-X", "PUT", "-H", "Content-Type: application/octet-stream", "--upload-file", str(tmp),
                                    "-w", "\n%{http_code}", f"{state['bucket']}/{name}"], capture_output=True, text=True, timeout=3600)
            finally:
                os.unlink(cfg.name)
            out = r.stdout.rsplit("\n", 1); code = out[1] if len(out) == 2 else "-1"
            try: srv = json.loads(out[0]).get("checksum", "").replace("md5:", "")
            except Exception: srv = ""
            if code in ("200", "201") and srv == md5:
                verdict = f"uploaded, md5 == (attempt {attempt})"; break
            print(f"{name}: attempt {attempt} HTTP {code} server md5 {srv[:8] or '-'} != local", file=sys.stderr, flush=True); time.sleep(30 * attempt)
        if verdict is None:
            if (have := server_files().get(name)) and have.get("checksum", "").replace("md5:", "") == md5: verdict = "landed despite the error, md5 =="
            else: print(f"{name}: FAILED after 5 attempts", file=sys.stderr, flush=True); sys.exit(1)
finally:
    tmp.unlink(missing_ok=True)
secs = time.time() - t0
with open(table, "a") as t: t.write(f"{name}\t{size}\t{sha}\t{md5}\t{md5}\t{secs:.0f}\n")
print(f"{name}\t{size}\tsha256 {sha[:16]}...\tmd5 {md5}\t{verdict}\t{secs:.0f} s", flush=True)
