# SPDX-License-Identifier: MIT
"""zenodo_upload.py <metadata.json> <files dir> <state.json> - draft Zenodo deposition (never publishes):
create the deposition with a prereserved DOI (or reuse the one in state.json), set metadata, upload the files of
zenodo/FILES.md through the bucket API one by one (streamed; up to 5 attempts each; a file whose md5 already matches on
the server is skipped), compare every reported checksum with the local md5. Token: ~/.zenodo_token (never printed).
Sandbox: ZENODO_BASE=https://sandbox.zenodo.org. Prints only ids, DOI, file names, sizes and checksum verdicts."""
import hashlib, json, os, sys, time
from pathlib import Path
import urllib.request, urllib.error

BASE = os.environ.get("ZENODO_BASE", "https://zenodo.org")
TOKEN = Path("~/.zenodo_token").expanduser().read_text().strip()
FILES = ["cohort558_q4k.tar", "cohort558_q4k.tar.sha256", "cohort558_q16k.tar", "cohort558_q16k.tar.sha256",
         "MANIFEST.tsv", "SHA256SUMS.q4k", "SHA256SUMS.q16k", "DATASET.md", "FORMAT.md"]


def api(method, url, data=None, ctype="application/json", stream=None, length=None):
    req = urllib.request.Request(url, method=method)
    req.add_header("Authorization", "Bearer " + TOKEN)
    if data is not None:
        body = json.dumps(data).encode(); req.add_header("Content-Type", ctype); req.data = body
    if stream is not None:
        req.add_header("Content-Type", "application/octet-stream"); req.add_header("Content-Length", str(length)); req.data = stream
    try:
        with urllib.request.urlopen(req, timeout=3600) as r:
            txt = r.read().decode()
            return r.status, (json.loads(txt) if txt else {})
    except urllib.error.HTTPError as e:
        txt = e.read().decode(errors="replace")
        return e.code, {"error": txt[:800]}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()


def main():
    meta_path, d, state_path = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    meta = json.load(open(meta_path))
    state = json.load(open(state_path)) if state_path.exists() else {}
    if "id" not in state:
        code, r = api("POST", f"{BASE}/api/deposit/depositions", {"metadata": {"prereserve_doi": True}})
        if code not in (200, 201): sys.exit(f"create failed: {code} {r}")
        state = {"id": r["id"], "doi": r["metadata"]["prereserve_doi"]["doi"], "bucket": r["links"]["bucket"], "files": {}}
        json.dump(state, open(state_path, "w"), indent=1)
    print("deposition", state["id"], "reserved DOI", state["doi"])
    code, r = api("PUT", f"{BASE}/api/deposit/depositions/{state['id']}", {"metadata": {**meta, "prereserve_doi": True}})
    if code != 200: sys.exit(f"metadata failed: {code} {r}")
    print("metadata set:", r["metadata"].get("title", "")[:60], "| doi", r["metadata"].get("prereserve_doi", {}).get("doi"))
    code, listed = api("GET", f"{BASE}/api/deposit/depositions/{state['id']}/files")
    remote = {f["filename"]: f for f in listed} if code == 200 and isinstance(listed, list) else {}
    for name in FILES:
        p = d / name; size = p.stat().st_size; local = md5(p)
        if name in remote and remote[name].get("checksum", "").replace("md5:", "") == local:
            print(f"{name}\t{size}\tmd5 {local}\talready on server, checksum ==")
            state["files"][name] = {"bytes": size, "md5": local, "server": "md5:" + local}; json.dump(state, open(state_path, "w"), indent=1); continue
        ok = False
        for attempt in range(1, 6):
            with open(p, "rb") as f:
                code, r = api("PUT", f"{state['bucket']}/{name}", stream=f, length=size)
            if code in (200, 201) and r.get("checksum", "").replace("md5:", "") == local and r.get("size") == size:
                print(f"{name}\t{size}\tmd5 {local}\tuploaded, checksum == (attempt {attempt})"); ok = True
                state["files"][name] = {"bytes": size, "md5": local, "server": r.get("checksum")}; json.dump(state, open(state_path, "w"), indent=1); break
            print(f"{name}\tattempt {attempt} failed: {code} {str(r)[:200]}"); time.sleep(30 * attempt)
        if not ok: sys.exit(f"upload failed: {name}")
    code, r = api("GET", f"{BASE}/api/deposit/depositions/{state['id']}")
    print("state:", r.get("state"), "submitted:", r.get("submitted"), "files on server:", len(r.get("files", [])))


if __name__ == "__main__":
    main()
