"""Refuse a stale TestPyPI upload: installed package bytes must match this run's wheel."""
from importlib.metadata import distribution
from pathlib import Path
import hashlib
import json
import sys
import zipfile

wheels = list(Path(sys.argv[1]).glob("*.whl"))
assert len(wheels) == 1, wheels
dist = distribution("panvram")
checked = 0
with zipfile.ZipFile(wheels[0]) as archive:
    for name in archive.namelist():
        if name.endswith("/") or ".dist-info/" in name:
            continue
        expected = archive.read(name)
        actual = Path(dist.locate_file(name)).read_bytes()
        assert actual == expected, f"installed wheel mismatch: {name}"
        checked += 1
assert checked > 0
print(json.dumps({"installed_wheel": "PASS", "files": checked, "version": dist.version,
                  "wheel_sha256": hashlib.sha256(wheels[0].read_bytes()).hexdigest()}))
