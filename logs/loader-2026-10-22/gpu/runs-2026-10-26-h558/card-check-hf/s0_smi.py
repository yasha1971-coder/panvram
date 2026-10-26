import subprocess, time, datetime
T_STAGE = {}
def stage(name, t0): T_STAGE[name] = time.time() - t0; print(f"STAGE {name}: {T_STAGE[name]:.1f} s  (utc {datetime.datetime.utcnow():%H:%M:%S})", flush=True)
def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True, executable="/bin/bash")
    if r.returncode: print(f"$ {c}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}", flush=True); raise RuntimeError(c)
    return r.stdout
print(sh("nvidia-smi --query-gpu=name,memory.total,clocks.max.sm,driver_version --format=csv").strip(), flush=True)
