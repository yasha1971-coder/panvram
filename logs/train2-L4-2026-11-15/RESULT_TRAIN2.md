# TRAIN-2 — RESULT (panvram cohort vs reference-only, stratified held-out eval)

Pre-registration: ACEAPEX/train1/PREREG_TRAIN2.md (738990d), committed BEFORE the run.
panvram FROZEN 01f0ce11 (code UNCHANGED). Model/hypers from examples/train_small.py:
TinyCNNLM d=128 L=8, W=2048, steps=2000, batch=64, lr=2e-3, holdout=50 (split by assembly),
reference CHM13 t2t.fa, seeds 20261115 / 1 / 2. Weights saved per (seed, model). GPU: L4
(self-deleting GCP VM), panvram built from pinned source 01f0ce1 (PANVRAM_CUDA=1).

## Eval sets (bits/base on the 50 held-out assemblies, never trained on)
- (a) random — `ev.sample` (TRAIN-1 replication / control).
- (b) HLA(chr6 MHC)+KIR(chr19)+SegDup — assembly windows whose CHM13 orthologue (release-2
  assembly-vs-CHM13 PAF `hprc465vschm13.aln.paf.gz`) falls in CHM13 chr6:28,000,000-34,000,000
  (MHC/HLA), chr19:54,000,000-56,000,000 (KIR/LRC), or `chm13v2.0_SD.bed` (T2T CHM13v2.0
  segmental duplications). 1024 windows/seed, fetched via panvram.Cohort.fetch(..., tokens=True).
  Method note: panvram fetch is in ASSEMBLY coordinates; the CHM13 loci were mapped to each
  held-out assembly's contigs through the PAF (stated, not a native CHM13-coord fetch).
- (c) SV-insertion — assembly windows with NO CHM13 alignment in the PAF (inserted vs CHM13).
  1024 windows/seed, same fetch.

## Prediction -> fact (mean over 3 seeds of reference_only - cohort, bits/base; >0 = cohort better)
| set | prediction (prereg) | fact mean | per-seed (default / 1 / 2) | verdict (prereg rule) |
|---|---|---|---|---|
| random         | ~0 (\|Δ\|<0.01) | **-0.0022** | -0.0012 / -0.0051 / -0.0002 | no effect |
| HLA+KIR+SegDup | cohort better ≥ +0.05 | **-0.0045** | -0.0001 / -0.0104 / -0.0030 | no effect (cohort slightly WORSE; prediction NOT met) |
| SV-insertion   | cohort better ≥ +0.05 | **+0.0302** | +0.0482 / +0.0232 / +0.0191 | weak/inconclusive (prediction NOT met) |

Absolute bits/base (cohort / reference_only), per set, mean of the three seeds:
random 1.721/1.719; HLA+KIR+SegDup 1.705/1.700; SV-insertion 1.685/1.715.

## Decision rule (verbatim from prereg 738990d, NOT reformulated)
"per set ≥ +0.05 → effect; 0.01 ≤ set < 0.05 → weak/inconclusive; < 0.01 → no effect. The
hypothesis is confirmed only if BOTH HLA+KIR+SegDup AND SV-insertion reach ≥ +0.05, while
random stays < 0.01."

Applying it: random |−0.0022| < 0.01 ✓ (control clean, as TRAIN-1 +0.0022). HLA+KIR+SegDup
−0.0045 → no effect (NOT ≥ +0.05). SV-insertion +0.0302 → weak/inconclusive (NOT ≥ +0.05).

## VERDICT: hypothesis NOT confirmed (hypothesis_confirmed = false).
The cohort does NOT reach the pre-registered +0.05 bit/base bar on either stratified set.
HONEST note on the one directional signal: on **SV-insertion windows the cohort is better on
all three seeds** (+0.048 / +0.023 / +0.019; mean +0.030), i.e. a consistent-sign but
sub-threshold effect — the seed-default +0.0482 sits right at the 0.05 bar but is not
reproduced by seeds 1/2, so by the fixed rule it is weak/inconclusive, not a win. HLA+KIR+
SegDup shows no benefit (slightly negative). This refines TRAIN-1's global null: the only
place the 558-assembly cohort trends toward helping is sequence INSERTED relative to CHM13,
and even there the gain is below the pre-registered significance bar at this scale
(W=2048, 2000 steps, 50-assembly holdout).

Logs/weights: ACEAPEX/train1/logs/train2/ (train2_summary.json, train2.seed-*.json,
weights/{cohort,reference_only}.seed-*.pt, driver.log, STATUS). Kit: pubrepo train1-kit
(train2_gcloud.sh / train2_driver.sh / train2_eval.py / train2_precompute_windows.py). VM self-deleted.
