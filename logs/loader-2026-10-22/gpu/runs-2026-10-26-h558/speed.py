import subprocess, time, hashlib, concurrent.futures as cf
URL = "https://zenodo.org/records/23232317/files/cohort558_q16k.tar.part00"
N = 524288000; K = 4; step = N // K
def get(i):
    a = i * step; b = N - 1 if i == K - 1 else a + step - 1
    subprocess.run(["curl", "-fsSL", "--retry", "5", "-r", f"{a}-{b}", "-o", f"/content/p.{i}", URL], check=True)
t = time.perf_counter()
with cf.ThreadPoolExecutor(K) as ex: list(ex.map(get, range(K)))
dt = time.perf_counter() - t
h = hashlib.sha256()
for i in range(K): h.update(open(f"/content/p.{i}", "rb").read())
ok = h.hexdigest() == "f7df5843b5eee7b88f0a34df86c10d8650ebe9f85ca8ee872da634a262fe9bf4"
print(f"part00 {N} B in {dt:.1f} s = {N/dt/1e6:.1f} MB/s, 4 range streams, sha256 ok {ok}; "
      f"16 parts at this rate: {16*dt/60:.1f} min (8 098 477 149 B: {8098477149/(N/dt)/60:.1f} min)", flush=True)
