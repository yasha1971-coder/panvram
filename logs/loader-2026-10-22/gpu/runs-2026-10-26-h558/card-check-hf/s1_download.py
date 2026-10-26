t0 = time.time()
B = "https://huggingface.co/datasets/Yakiv-Contour/panvram-hprc558/resolve/main"
sh("mkdir -p /content/dl /content/cohort")
sh(f"cd /content/dl && curl -fsSL -o PARTS.sha256 {B}/PARTS.sha256 && curl -fsSL -o tar.sha256 {B}/cohort558_q16k.tar.sha256 && "
   f"seq -w 0 15 | xargs -P 4 -I{{}} curl -fsSL --retry 5 -o cohort558_q16k.tar.part{{}} {B}/cohort558_q16k.tar.part{{}}")
print(sh("du -sb /content/dl").strip(), flush=True)
stage("download_hf_16_parts_4_streams", t0)
