# TRAIN-1 — RESULT (panvram cohort vs reference-only, holdout bits/base)

Date 2026-11-15. GPU: NVIDIA L4 (24 GB), self-deleting GCP VM (g2-standard-16, <zone>),
DL image pytorch-2-9-cu129-ubuntu-2204-nvidia-580. panvram built from pinned source
01f0ce11 (PANVRAM_CUDA=1), code UNCHANGED. Config: examples/train_small.py, W=2048,
steps=2000, batch=64, holdout=50, cohort 508 train / 50 heldout assemblies, reference=CHM13 t2t.fa.
Pre-registration: ACEAPEX/train1/PREREG_TRAIN1.md (1e9505b), committed before the run.

| seed | cohort bit/base | reference-only bit/base | cohort − ref | verdict |
|---|---|---|---|---|
| default | 1.7248 | 1.7236 | +0.0012 (cohort worse) | no effect |
| 1 | 1.7314 | 1.7263 | +0.0051 (cohort worse) | no effect |
| 2 | 1.7082 | 1.7080 | +0.0002 (cohort worse) | no effect |

Mean cohort−ref = +0.0022 bit/base (cohort slightly HIGHER = worse on every seed).
Pre-registered rule: cohort must be >= 0.01 bit/base LOWER on holdout to count as an effect;
< 0.01 => no effect.

**VERDICT: NO EFFECT.** At this scale (W=2048, 2000 steps, 50-assembly holdout) the 558-assembly
cohort does NOT improve held-out bits/base over reference-only; it is marginally worse on all
three seeds, well inside the 0.01 bit/base no-effect band. Honest negative result — the first
real TRAIN-1 number. tok/s ~1.13M (L4; GPU type affects only speed, not the bits/base verdict).

Logs: ACEAPEX/train1/logs/ (train1_summary.json, train1_runs.jsonl, STATUS, driver.log).
Kit: pubrepo branch train1-kit (train1_gcloud.sh, train1_driver.sh); VM self-deleted.
