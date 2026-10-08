#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""check_source_links.py <aceapex clone> [<ref>] - every link of README.md, EVIDENCE.md and DATASET.md of the form
https://github.com/yasha1971-coder/aceapex/tree/<tag>/<path> must point at a path that exists in <ref> (default: the tag
named in the link) of the aceapex clone (`git cat-file -e <ref>:<path>`). Exit 1 if any link is broken."""
import re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
clone = sys.argv[1]; ref_override = sys.argv[2] if len(sys.argv) > 2 else None
pat = re.compile(r"https://github\.com/yasha1971-coder/aceapex/tree/([A-Za-z0-9._-]+)(?:/([^\s)`>\"|;,]+))?")
seen, broken = {}, []
for name in ("README.md", "EVIDENCE.md", "DATASET.md"):
    for ln, line in enumerate((ROOT / name).read_text().splitlines(), 1):
        for tag, path in pat.findall(line):
            if path in ("", "<path>"): continue      # bare tag link / placeholder in the explanation
            key = (ref_override or tag, path)
            if key not in seen:
                seen[key] = subprocess.run(["git", "-C", clone, "cat-file", "-e", f"{key[0]}:{key[1]}"], capture_output=True).returncode == 0
            if not seen[key]: broken.append(f"{name}:{ln} {key[0]}:{key[1]}")
print(f"links checked: {len(seen)} distinct paths; broken: {len(broken)}")
for b in broken: print("BROKEN", b)
sys.exit(1 if broken else 0)
