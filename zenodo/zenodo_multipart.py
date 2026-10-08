# SPDX-License-Identifier: MIT
"""zenodo_multipart.py <record id> <file> [part MiB] - upload one large file to a Zenodo DRAFT through the InvenioRDM
multipart transfer (type "M"): initiate, PUT every part with curl (token through a 0600 config file, 5 attempts per
part), commit, then compare the server checksum with the local md5. Never publishes. Fallback for files that the
legacy single-PUT bucket API answers with 502."""
import hashlib, json, os, subprocess, sys, tempfile, time
from pathlib import Path
import urllib.request, urllib.error

BASE = os.environ.get("ZENODO_BASE", "https://zenodo.org")
TOKEN = Path("~/.zenodo_token").expanduser().read_text().strip()


def api(method, url, data=None):
    req = urllib.request.Request(url, method=method); req.add_header("Authorization", "Bearer " + TOKEN)
    if data is not None: req.add_header("Content-Type", "application/json"); req.data = json.dumps(data).encode()
    try:
        with urllib.request.urlopen(req, timeout=600) as r:
            t = r.read().decode(); return r.status, (json.loads(t) if t else {})
    except urllib.error.HTTPError as e: return e.code, {"error": e.read().decode(errors="replace")[:600]}


def curl_put_range(url, path, offset, length):
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(Path(path).parent), prefix=".curlcfg", suffix=".tmp") as cfg:
        os.chmod(cfg.name, 0o600); cfg.write(f'header = "Authorization: Bearer {TOKEN}"\n')
    try:
        # stream exactly one part: dd from the file into curl's stdin
        dd = subprocess.Popen(["dd", f"if={path}", "bs=1M", f"skip={offset // (1 << 20)}", f"count={(length + (1 << 20) - 1) // (1 << 20)}", "status=none"], stdout=subprocess.PIPE)
        head = subprocess.Popen(["head", "-c", str(length)], stdin=dd.stdout, stdout=subprocess.PIPE); dd.stdout.close()
        r = subprocess.run(["curl", "-sS", "-K", cfg.name, "-X", "PUT", "-H", "Content-Type: application/octet-stream", "-H", f"Content-Length: {length}",
                            "--data-binary", "@-", "-w", "\n%{http_code}", url], stdin=head.stdout, capture_output=True, text=True, timeout=4 * 3600)
        head.stdout.close(); dd.wait(); head.wait()
    finally:
        os.unlink(cfg.name)
    out = r.stdout.rsplit("\n", 1); return (int(out[1]) if len(out) == 2 and out[1].isdigit() else -1), (out[0][:300] if out else r.stderr[:300])


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()


def main():
    rid, path = sys.argv[1], Path(sys.argv[2]); part = int(sys.argv[3]) << 20 if len(sys.argv) > 3 else 512 << 20
    size = path.stat().st_size; parts = (size + part - 1) // part; key = path.name
    files_url = f"{BASE}/api/records/{rid}/draft/files"
    code, d = api("GET", files_url); entries = {e["key"]: e for e in d.get("entries", [])} if code == 200 else {}
    if key in entries and entries[key].get("status") == "completed":
        print(key, "already completed on the server:", entries[key].get("checksum")); return
    if key in entries:
        api("DELETE", f"{files_url}/{key}")
    code, d = api("POST", files_url, [{"key": key, "size": size, "transfer": {"type": "M", "parts": parts, "part_size": part}}])
    if code not in (200, 201): sys.exit(f"initiate failed: {code} {d}")
    entry = [e for e in d["entries"] if e["key"] == key][0]; links = entry["links"]
    part_links = links.get("parts")
    if not part_links: sys.exit(f"no multipart links returned (transfer type M not available): {json.dumps(entry)[:400]}")
    print(f"{key}: {size} B in {parts} parts of {part >> 20} MiB; multipart initiated", flush=True)
    for p in part_links:
        n = p["part"]; off = (n - 1) * part; ln = min(part, size - off)
        for attempt in range(1, 6):
            c, body = curl_put_range(p["url"], path, off, ln)
            if c in (200, 201, 204): print(f"  part {n}/{parts} ok ({ln} B, attempt {attempt})", flush=True); break
            print(f"  part {n}/{parts} attempt {attempt} failed: {c} {body[:120]}", flush=True); time.sleep(30 * attempt)
        else: sys.exit(f"part {n} failed")
    code, d = api("POST", links["commit"])
    if code not in (200, 201): sys.exit(f"commit failed: {code} {d}")
    local = md5(path); server = d.get("checksum", "")
    print(key, "committed; server", server, "| local md5:" + local, "| ==", server == "md5:" + local, flush=True)
    if server != "md5:" + local: sys.exit("checksum mismatch")


if __name__ == "__main__":
    main()
