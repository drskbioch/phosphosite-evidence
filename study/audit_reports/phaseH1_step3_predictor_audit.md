# Phase H1, Step 3 — Predictor Audit + Gate G3

**Scope:** PART 1 ONLY (Task 0 + Task A) of Step 3. This is the tier-deciding step (Tier A vs. Tier B, per H11) — mandatory stop after Task A for explicit sign-off before any of Task B–E runs. Not Step 4, not the manuscript, not a revival of the temporal holdout (G2 failed decisively in Step 2 — see `phaseH1_step2_temporal.md`).
**Date:** 2026-09-25

---

## Task 0 — Setup and state check

- `study/processed_data/epsd/sei/sei_protein.parquet` confirmed present and readable (2,170,157 bytes).
- `study/processed_data/epsd/temporal/site_first_year.parquet` confirmed present and readable (15,911,106 bytes).
- **G1/G2 recorded outcomes, quoted from the source reports (re-read in full before starting this brief, per standing rule 1):**
  - G1 (`phaseH1_sei_matching.md`, Task D): *"G1: PASSES. Disorder survives matching. Current framing holds, strengthened — not merely holds: the matched effect (0.042) is larger than the raw unmatched pooled baseline (0.031)."*
  - G2 (`phaseH1_step2_temporal.md`, Task D): *"G2: FAILS. Not usable... Finding 6 is extended to 'both standard partitionings fail'... Step 3 uses SEI stratification only, no temporal component."*
- `git status`: same two benign items as Steps 1–2 (`frontend` gitlink, `legacy/deploy/deepphos-frontend` embedded-repo content) plus the harness lock file. No new dirty state. Not stopping, per the same established reasoning.
- Created `study/scripts/earlier_predictor_audit/step3/` and `study/processed_data/epsd/predictor_audit/`.

---

## Task A — Candidate survey and design risks

### A.1 — Tool candidates against H7's rule

All entries below verified via web search this session (not recalled from training-data memory alone), sources listed per tool. **Not pre-filtered** — reporting what's actually available, H7's rule applied after.

| Tool | Weights downloadable? | CPU-runnable at 5-10k sites? | Batch server alternative | License/terms for reporting results | Species support |
|---|---|---|---|---|---|
| **MusiteDeep** | **Y** — GitHub (`duolinwang/MusiteDeep`), pretrained models included in repo | **Likely Y**, but the original repo targets Python 2.7 + Keras 1.1.0/Theano — a legacy stack; an updated repo (`duolinwang/MusiteDeep_web`) exists and should be checked for a modern-dependency version before committing to this tool. Not yet verified hands-on. | Yes, a webserver also exists (published in *NAR* 2020, `W1`) | Open, academic publication (Bioinformatics 2017 + NAR 2020); no restrictive license found | **Human-only** — repo's own description: "At present, MusiteDeep only provides prediction of human phosphorylation sites" |
| **GPS 6.0** | **N** — server-only (`gps.biocuckoo.cn`), no bulk weights-download found | N/A (not locally runnable) | **Y** — web server with multiple prediction modules (LightGBM fast mode, integrated DNN mode); batch-size limits not confirmed, needs a direct check before relying on it | Academic publication (*NAR* 2023); server ToS for bulk/scripted queries not confirmed — needs checking before use | Not species-restricted in what was found; trained on a large pooled p-site set, kinase-specific models span many kinases |
| **NetPhos 3.1 / 4.1** | **Possibly Y** — DTU Health Tech lists an "Academic Download" page for NetPhos 3.1 (Linux package), distinct from the public web server | Not verified hands-on; DTU academic-download packages are typically standalone C/Perl binaries, historically CPU-light | Y — public web server (`services.healthtech.dtu.dk`), per-submission size limits not confirmed | **Academic-use license required** (DTU's standard academic-download terms) — needs explicit check on whether benchmark-result reporting is permitted under those terms before use | General eukaryotic, not species-restricted |
| **PhosphoSVM** | **Y** — "server, trained model, and all datasets" stated available at `sysbio.unl.edu/PhosphoSVM` | Likely Y (SVM + sequence-attribute features, not deep learning — should be CPU-cheap), not verified hands-on | Also has a web server | Academic publication (*Amino Acids* 2014); license terms on the host site not confirmed | Not species-restricted in what was found |
| **DeepPhos (Luo et al. 2019)** | **Y** — GitHub (`USTC-HIlab/DeepPhos`, published in *Bioinformatics* 2019) | Not verified hands-on; densely-connected CNN, plausibly CPU-feasible at this sample size but not confirmed | Not found | Open, academic publication; **name collision with this project's own former legacy branding** — H14(b) already requires removing "DeepPhos" naming from this project regardless of whether this tool is used | Not confirmed multi-species; PSP-derived training data (typically human-dominant) |
| **Musite (original, 2010)** | **Y** — `musite.sourceforge.net`, **GNU GPL**, Java, source + models | Likely Y (2010-era classifier, lightweight relative to modern deep nets) | Not found (this is the offline tool, predates the modern web-server-first pattern) | **GPL — clean for reporting results** | **Multi-species by design** — "general prediction models were trained for six organisms" |
| **PLM-based (ESM-2/ProtBERT head)** | Mixed — no single canonical published weights found; multiple 2024–2025 preprints (PhosF3C, a kinase-specific ESM-2-decoder model) exist but none confirmed to ship ready-to-run inference weights as of this search | Depends on which; ESM-2 650M itself is the model already used for this project's own legacy PhosConsensus-Predict baseline (H7 names this as a baseline, not a contender) | N/A | Varies by paper | Not applicable / varies |
| **DeepPSP, PhosIDN, Quokka** (additional established tools found during the search) | DeepPSP/PhosIDN: GitHub repos exist (`ustchangyuanyang/PhosIDN` for PhosIDN); Quokka: logistic-regression-based, smaller/older, download status not confirmed | Not verified hands-on for any of the three | Not confirmed | Academic publications; license terms not confirmed | Not confirmed for any of the three |

**Honest gaps in this table, stated plainly rather than papered over:** several "Y" entries for downloadability are based on the tool's GitHub/site *existing* and being described as containing weights/models, not on this session actually cloning and running each one — that hands-on verification is Task B/C work, gated on sign-off. Several license/ToS questions (GPS 6.0 and NetPhos server scripted-query terms, specifically) are flagged as **not yet confirmed** and need resolving before those tools are used, not assumed permissive.

### A.2 — Temporal confound check on the stratification covariate

Computed directly (`site_first_year.parquet` × `sei_protein.parquet`, n=1,439,758 sites with both values, mapped proteins only):

- **Spearman ρ = -0.0112** (p=2.1e-41), **Pearson r = -0.0305** (p=8.7e-294).
- Statistically significant only because of the huge sample size — **the effect size is negligible**, essentially zero correlation between SEI_primary and first-reported year at the site level.
- Within-tier: GOLD ρ=-0.077, SILVER ρ=-0.064, BRONZE ρ=+0.045 — all small.

**Interpretation, not a correction attempt (per this task's own instruction):** SEI_primary is **not** meaningfully era-correlated, unlike the tier variable itself (Step 2's ~14-year GOLD/BRONZE gap). This means an SEI-decile-balanced evaluation sample will **not** be implicitly era-stratified in any material way — a predictor-gap-by-SEI-decile result in Task C is unlikely to be confounded by "predictors trained on more recent data performing differently on newer sites," since SEI-decile membership doesn't track era. This is reassuring for interpreting Task C's results cleanly, though the correlation being non-zero and statistically significant (given n) means it's not literally zero — reported precisely, not rounded to "no relationship."

### A.3 — Negative-sampling risk (reported and flagged, not resolved)

**The risk, stated plainly:** negatives are non-annotated S/T/Y positions on SEI-matched proteins, explicitly `label_type: unlabelled` per H6 — never asserted as true negatives. On **understudied (low-SEI) proteins**, a higher fraction of these "unlabelled" positions are plausibly real, undiscovered phosphosites — simply because fewer of that protein's true sites have been found at all (fewer studies → fewer sites detected → more of the protein's real phosphosites remain in the "unlabelled" pool, undetected, by construction). On **well-studied (high-SEI) proteins**, more of the true phosphosite landscape has already been found and moved into the positive (GOLD/BRONZE) set, so the "unlabelled" pool there is more likely to be genuinely non-phosphorylated.

**Why this matters for G3's interpretation, in both directions:** if a predictor's AUROC looks *worse* on low-SEI proteins, that could genuinely mean the predictor performs worse there (real signal) — **or** it could mean the negative set is more contaminated there with undiscovered true positives, which would *deflate* the AUROC by penalizing the predictor for correctly scoring likely-real sites high (while the "unlabelled" label scores them as wrong). **This cuts in the opposite direction from a naive "predictors underperform on understudied proteins" reading** — contamination would make predictors look *artificially worse* specifically where the paper's thesis expects a real gap, which could inflate an apparent understudied-vs-well-studied AUROC gap that isn't (only) about predictor behavior.

**Not resolved computationally, per this task's explicit scope.** Reported here so it can travel as a mandatory caveat with whatever G3 finds, per the brief's own Task D requirement — regardless of which way the result comes out.

### A.4 — Scope recommendation

**Recommendation: human-only**, driven directly by A.1's findings — MusiteDeep, the single clearest H7-qualifying tool found (downloadable weights, GitHub, established, clean license), is **explicitly human-only** by its own documentation. A mammalian-pooled scope would either exclude MusiteDeep entirely (weakening the shortlist against H7's "≥4 external tools" target) or force a species-inconsistent comparison (some tools scored only on human sites, others on all three). None of the other surveyed tools were confirmed to have clean, verified multi-species support at the downloadable-weights level except the original 2010 Musite (six organisms, but a much older/weaker baseline architecture) — the modern, strongest H7-qualifying candidates skew human-focused.

**Axis 1's human-only limitation (Step 1) was for a different reason (raw annotation files gone) and is unrelated here** — most of these predictors are sequence/PLM-based and don't consume domain annotations directly, so Axis 1's specific data-availability gap doesn't propagate to this step. The human-only recommendation here is driven by *tool* species-support, not by any shared data limitation with Step 1.

---

## Task A — STOP. Sign-off needed on four items before Task B

1. **Final tool list.** Candidate shortlist meeting H7 pending hands-on verification: **MusiteDeep, PhosphoSVM, Musite (original), DeepPhos (Luo et al.)** as the four most likely to clear H7 cleanly (downloadable weights, reasonable license clarity) — GPS 6.0 and NetPhos likely need the batch-server path instead of local weights, pending their ToS/limit checks. Confirm this shortlist, or adjust it.
2. **Human-only vs. mammalian-pooled scope.** Recommendation above: human-only, driven by MusiteDeep's species restriction. Confirm or override.
3. **Whether A.2's temporal-correlation finding changes the stratification design.** Recommendation: no change needed — the correlation is negligible, SEI-decile stratification is not implicitly era-stratified. Confirm.
4. **A.3's negative-sampling caveat — carry forward as-is, or attempt a mitigation?** No mitigation was invented (per instruction not to invent a signal that doesn't exist in the data) — if a genuine "minimum evidence of absence" signal exists (e.g., a site with multiple independent MS runs reporting *no* modification at that position, if EPSD or its sources capture that), it would need to be identified in the data before being proposed, not assumed. Confirm the caveat travels as-is into Task D, or provide direction on what to check for a possible mitigation.

---

## Follow-up round (2026-09-25) — MusiteDeep feasibility, license split, backup candidate, A.4 recheck

### 1. MusiteDeep environment feasibility — tested directly, not estimated

**Checked this machine directly:** no `conda`, no Python 2.7, no Docker installed (`which`/`where.exe` all negative). Only Python 3.12/3.13 available via the `py` launcher.

**Fetched both MusiteDeep variants' actual requirements** (not assumed from the earlier search snippet):
- Main repo: Python 2.7, Keras 1.1.0, Theano. Docs state *"tested in Linux and Mac OS X"* — **no Windows mention at all.**
- `MusiteDeep_Keras2.0` subfolder (the one that looked like a possible Python 3 escape route): **still Python 2.7** — it only upgrades Keras (2.1.2) and TensorFlow (1.3.0) on top of the same Python 2.7 base, confirmed by fetching that subfolder's own README directly. Not a Python 3 path at all.
- Pretrained weights are **not bundled** in either variant — a separate download step, destination unconfirmed.

**Verdict: this exceeds the 1-2 hour threshold — saying so now, as asked.** TensorFlow 1.3.0 + Python 2.7 on Windows, with no conda/Docker pre-installed and no official Windows support claimed by the tool itself, is a real, multi-hour-or-worse risk (likely needs Docker Desktop installed first — itself a nontrivial system change requiring admin/restart — plus locating or building a compatible legacy image, since modern pip cannot reliably resolve TF 1.3.0 wheels for Python 2.7 in 2026). **Recommend deprioritizing MusiteDeep as the first tool attempted in Task B**, or explicitly budgeting real setup time for it with the backup (item 3) ready if it stalls.

### 2. License split — (a) weight redistribution vs. (b) run-and-report, per tool

**Clarifying context that changes the framing:** H10 only deposits *this project's own* package/SEI table/evaluation sample/corpus on Zenodo — it does not call for redistributing any external tool's weights. So **(a) may not even be operationally needed for any of these tools** — worth confirming, but if true, (a)'s uncertainty is moot for every tool below except as a matter of completeness.

| Tool | (a) Weight redistribution | (b) Run locally + report results |
|---|---|---|
| MusiteDeep | GPL-2.0 — **permitted**, propagation terms apply if redistributed | **Standard academic practice, no issue** |
| Musite (original) | GPL — **permitted**, same as above | **No issue** |
| PhosphoSVM | Not confirmed on the host site — genuinely open | **No issue** — running public code and reporting its output in a benchmark table doesn't require explicit permission under normal academic norms, independent of (a)'s status |
| DeepPhos | **No LICENSE file found at all** in the repo — default copyright applies, redistribution genuinely unclear without contacting the authors | **No issue**, same reasoning as above — (b) doesn't depend on (a) being resolved |
| Phosformer (backup) | LICENSE file exists (exact type not confirmed from the fetch) | **No issue** |

**Bottom line: (b) is not actually blocked for any tool on this list.** The only real open question is (a), and only for DeepPhos is it genuinely unresolved (no license at all) rather than just unconfirmed-but-likely-fine — and per the H10 clarification above, (a) likely isn't needed regardless.

### 3. Backup candidate: **Phosformer** (`esbgkannan/phosformer`, *Bioinformatics* 2023)

Modern, PyTorch, auto-downloads weights from HuggingFace Hub on first run (no manual weight-hunting), explicitly CPU-feasible per its own docs ("PyTorch should always work on CPU"). Lowest setup-friction candidate found this round.

**Caveat, stated honestly, not oversold:** Phosformer does **kinase-specific** scoring (paired kinase-domain + substrate-peptide input), not a direct general phosphosite/non-phosphosite binary classifier like the other four. It answers a related but different question. Using it as a drop-in for the general-predictor comparison would need a specific design choice in Task B (e.g., max score across a candidate kinase panel as a general "phosphorylatable" proxy) — not plug-and-play equivalent to MusiteDeep/PhosphoSVM/DeepPhos/Musite. Valuable specifically as a low-friction fallback if a primary tool stalls, not a silent equivalent substitute.

### A.4 recheck — is human-only really the tightest constraint, or just the simplest to state?

**It was entirely MusiteDeep-driven, and that changes now that MusiteDeep is a feasibility risk:**

| Tool | Species support, verified |
|---|---|
| PhosphoSVM | **Multi-species natively** — training data (Phospho.ELM-derived) explicitly includes *Homo sapiens* and *Mus musculus* |
| Musite (original) | **Multi-species by design** — general models trained for six organisms |
| DeepPhos | **Not documented** — README states no species scope at all; treat as unconfirmed, not verified either way |
| Phosformer (backup) | Human kinases only (its reference kinase-domain data is human-specific) |

**Revised picture:** if MusiteDeep is deprioritized per item 1, the remaining shortlist skews *toward* mammalian-pooled feasibility, not away from it — two of the other three tools already support mouse natively. Human-only is no longer the obviously tightest binding constraint; it was the simplest thing to state when MusiteDeep was still assumed primary. **This needs your call now, informed by this correction** — human-only remains the safer, more conservative choice (keeps DeepPhos's unconfirmed scope from being a problem, avoids any per-tool species inconsistency), but mammalian-pooled is more achievable than originally stated if that's preferred.

---

## Modern PLM-based general predictor search (2026-09-25)

Per sign-off requirement: search for ≥1 2021–2025 protein-language-model-based **general** (not kinase-specific) phosphosite predictor with downloadable weights, CPU/Kaggle-runnable, verified against sources not memory, reported before Task B starts.

Four candidates identified via web search and individually verified (GitHub repo and/or paper fetched directly, not taken from search-snippet summaries alone). One additional name (PhosphoLingo, bioRxiv 2022) surfaced but was not verified — flagged, not used. Two names from the earlier round were screened and excluded up front: **GPS-pPLM** (prokaryotic-specific — off-target for this human/mammalian project) and **NphosNet** (2026, confirmed via its own title to be N-phosphorylation-specific, a different PTM chemistry from the S/T/Y phosphorylation this project studies — off-target, excluded).

| Tool | General or kinase-specific | Weights downloadable | CPU feasibility | Species (training data) | License | Publication date (for leakage-free subset) |
|---|---|---|---|---|---|---|
| **PhosF3C** (`SkywalkerLuke/PhosF3C`) | **General S/T/Y** — confirmed, not kinase-specific (paper explicitly declines to train per-kinase-family models) | **Y** — paper (PMC12107248) states code+data "publicly accessible" at the GitHub repo; README requires downloading LoRA + Conformer checkpoint weights from a linked Google Drive into `model_ckp/`. Link presence confirmed; whether the Drive folder is currently populated was **not independently confirmed by downloading it** — needs a direct check in Task B before relying on it. | Config exposes `device: cuda` **or** `cpu` explicitly. ESM2-650M backbone with LoRA adapters — same embedding family already used by this project's own legacy PhosConsensus-Predict baseline, so local CPU/Kaggle feasibility is plausible and somewhat de-risked by prior project experience, but not benchmarked here. | **Multi-species**, not human-only — training data is "a comprehensive fusion of data sourced from MusiteDeep and DeepPSP," itself integrating UniProt/SwissProt, PhosphoSitePlus, and Phospho.ELM. Under a human-only scope decision, this is a training/eval population mismatch to disclose, same class of caveat as several existing shortlist tools. | **MIT** — confirmed via repo footer. Clean on both (a) and (b). | *Briefings in Bioinformatics* 26(3), bbaf242 — received 2025-02-14, accepted 2025-05-01, published May 2025; bioRxiv preprint 2024-12-25. **No explicit training-data collection cutoff stated in the paper** — recommend using the bioRxiv preprint date (2024-12-25) as the conservative publication-year proxy for the leakage-free-subset test, flagged as an approximation since the underlying MusiteDeep/DeepPSP/PSP/Phospho.ELM source data predates the paper itself. |
| **PTransIPs** (`StatXzy7/PTransIPs`) | **General S/T/Y** — separate Y and S/T models, no kinase-specific framing found | **Y**, but hosted on OneDrive (not GitHub releases/Zenodo) — less durable hosting than PhosF3C | `--device` defaults to `cuda:0` in the training script; no CPU feasibility statement found. Architecture (CNN+Transformer on ProtTrans+EMBER2 embeddings) is smaller than ESM2-650M, so CPU inference is plausible but **unverified**. | Not clearly documented; paper references SARS-CoV-2 phosphosites as a use case, general training scope otherwise undocumented in the fetched README | **No LICENSE file found** — same "run-and-report-only, no redistribution" treatment already applied to DeepPhos | arXiv 2308.05115, 2023-08-09 |
| **LMPhosSite** (`KCLabMTU/LMPhosSite`) | **General S/T/Y** | **Y** — Google Drive link, not in-repo, not GitHub releases | Trained on "200GB RAM, 18 CPUs, 1 Nvidia Tesla V100 16GB GPU," README recommends "similar configuration" — GPU implied preferred, ProtT5 embedding extraction is the heavy step; CPU feasibility unverified, likely slow but not confirmed infeasible. Older stack (Keras 2.4.3 / TF 2.3.1 / Python 3.8.3). | Not specified in README | **Apache-2.0** — confirmed LICENSE file present. Clean on both (a) and (b), cleanest redistribution story of the three PLM candidates. | *J. Proteome Research*, 2023-08-04 |
| **SAGEPhos** (`ZhangJJ26/SAGEPhos`) | **DISQUALIFIED — kinase-specific.** "Bio-Coupled Modal Fusion" explicitly distills "kinase sequence information" to refine substrate features — same category problem already ruled out for Phosformer (H7/A.4 backup candidate), for the same reason: a different task, adapting it would be an indefensible design choice. | Y (GitHub releases, `checkpoint.tar.gz`) — noted for completeness only | GPU implied strongly (2× A40 default in docs) | Not stated | Not stated (no LICENSE file found) | arXiv 2502.07384, 2025-02 (ICLR 2025) |
| PhosphoLingo (bioRxiv 2022) | Reportedly general, PLM-based | **Not verified this round** | **Not verified** | **Not verified** | **Not verified** | Not verified — flagged as a possible fifth option, not checked to the same depth as the four above due to time-boxing; can be verified in Task B if the shortlist needs a fifth tool |

**Recommendation, not a unilateral decision:** **PhosF3C** is the strongest candidate — explicitly general (not kinase-specific), cleanest license (MIT, permissive on both weight redistribution and run-and-report), most recent and directly relevant publication venue (*Briefings in Bioinformatics*, the same journal family this project targets), and architecturally continuous with this project's existing ESM2 experience from the legacy baseline. Its one open item is confirming the Google Drive checkpoint link is actually populated — a five-minute check to do at the start of Task B, not a blocker to reporting now. **LMPhosSite is the safer fallback** if PhosF3C's weights turn out to be unavailable or CPU inference proves infeasible — cleanest license (Apache-2.0, weights fully redistributable) though GPU-leaning. **PTransIPs is the weakest of the three viable options** — no license, less durable weight hosting, and no CPU feasibility evidence — usable under run-and-report-only terms if needed but not the first choice.

Per instruction, **not adding any of these to the tool list without your confirmation.** Awaiting your pick (PhosF3C recommended) before folding it into Task B's design, alongside the still-open Task A sign-off items above if not already resolved.

---

## Final verification round (2026-09-25) — PhosF3C checkpoint + MusiteDeep timeboxed attempt

### PhosF3C checkpoint — VERIFIED LIVE

Fetched the Google Drive folder directly (`https://drive.google.com/drive/folders/1I7qZw6cGUcVHmPYkt2hZ-tSsSiH-Gjlu`): loads without sign-in, shows two populated subfolders — `conformer` and `lora`, both modified **2024-12-17**, matching the README's expected `model_ckp/{lora,conformer}` structure and consistent with the paper's 2024-12-25 bioRxiv date. Not dead, not gated. **PhosF3C is locked in.** (Full end-to-end download-and-run was not performed here — that's Task B execution — but access-status verification, which is what was asked, is clean.)

### MusiteDeep — Kaggle timeboxed attempt: NOT EXECUTED, tooling gap on my end (not a MusiteDeep failure)

**Reported plainly, distinguishing this from a technical failure of the tool itself:** this session has no Kaggle CLI, no Kaggle API credentials (`~/.kaggle` absent, `KAGGLE_USERNAME`/`KAGGLE_KEY` unset), and no Kaggle-notebook-execution tool available to me. I cannot actually launch and run a Kaggle notebook from here — this is a capability gap in what I have access to in this environment, distinct from "MusiteDeep didn't produce predictions within 3 hours," which is the failure mode the standing instruction anticipated and which I have not observed because the attempt was never run. Flagging this explicitly rather than quietly reporting a pass/fail that didn't happen.

**The pre-authorized alternative path was checked instead, and it passed cleanly.** Fetched the MusiteDeep web server documentation directly from its primary source (*Nucleic Acids Research* 48(W1):W140–W146, 2020, `academic.oup.com/nar/article/48/W1/W140/5824154`), not a search-snippet summary:
- FASTA file upload up to **10 MB** per job
- Up to **5 concurrent jobs** per user
- Results retained 72 hours, up to 100MB per user
- Stated throughput: **"less than three minutes to predict for 1000 sequences per PTM type"**

A 10MB FASTA of typical peptide-length sequences comfortably holds several thousand sequences, and the stated ~3 min/1000-sequence throughput scales linearly to a fully usable turnaround for a 5,000–10,000-site evaluation sample (split across a handful of jobs if needed to stay under the 10MB/job cap). **This clears the "batch submission of a few thousand sequences" bar cleanly.**

**Recommendation:** use the web server as MusiteDeep's execution path for Task B, and treat the local Python-2.7/Theano Kaggle attempt as unnecessary rather than failed — same real MusiteDeep predictions, none of the legacy-stack risk already documented in the earlier feasibility check. If you'd rather have the Kaggle attempt actually run (e.g. to compare against the web server's output, or because the web server's specific PTM-type coverage differs from the local model), I can prepare a ready-to-paste Kaggle notebook script (Python 2.7 conda env creation, MusiteDeep clone, weight download, inference) for you to execute and report back — I just can't run it myself from this session.

### Finalized tool list (5 tools)

| Tool | Publication year (for leakage-free subset) | Execution path |
|---|---|---|
| PhosphoSVM | 2014 (*Amino Acids*) | Local / server |
| Musite (original) | 2010 | Local |
| DeepPhos (Luo et al.) | 2019 (*Bioinformatics*) — **run-and-report only, no weight redistribution** | Local |
| PhosF3C | 2024-12-25 (bioRxiv preprint date, used as conservative proxy — no explicit training-cutoff stated in the paper) | Local (CPU/GPU-configurable) |
| MusiteDeep | 2017 (*Bioinformatics*, original model) — note: web server itself published 2020 (*NAR*); using the 2017 model-publication date as the conservative (earlier) cutoff for leakage-free-subset purposes, since that's when the underlying model weights were trained | **Web server** (batch FASTA submission), not local, not Kaggle |

PTransIPs and LMPhosSite remain documented backups only, not in the initial run, per your instruction. Human-only scope confirmed. A.2/A.3 sign-offs stand unchanged.

**Both outcomes reported. Proceeding straight into Task B — no further gate, per instruction.**

---

## Task C, hands-on tool verification round (2026-09-25)

Per instruction, reporting after each tool clears or fails rather than batching. Task B's frozen sample is built (10,000 sites, checksum `0285523...477a6cb`); this section verifies each tool is *actually invocable*, not just downloadable, before running it against that sample.

### PhosphoSVM — **FAILS, documented tool failure**

Downloaded both files from `sysbio.unl.edu/PhosphoSVM/download.php` directly (site is HTTP-only, no TLS — required fetching raw HTML/files via curl rather than the standard fetch tool, which auto-upgrades to HTTPS and fails against this host):
- `svm-model.zip` (18.7MB): **raw libsvm model files only** (`S-model.txt`, `T-model.txt`, `Y-model.txt` — standard libsvm text format, RBF kernel, support vectors as ~36+ numbered feature indices). No scoring script, no feature-extraction code.
- `dataset.zip` (12.4MB): raw per-protein training sequences (P.ELM/PPA), not code either.
- The web prediction form (`prediction.php`) is **CAPTCHA-gated on every submission** ("input the correct four-letter word shown in figure") — not scriptable, and not viable for a 10,000-site batch regardless of automation (a CAPTCHA per submission is a manual bottleneck no matter how the sequences are batched).

**To actually score sequences with the downloaded model would require reverse-engineering PhosphoSVM's exact feature-encoding scheme from the paper's text and rebuilding that pipeline from scratch against someone else's trained weights.** Not attempted — a mismatched feature mapping would silently produce plausible-looking wrong scores, which is a worse outcome than an honest failure. **Recorded as a documented tool failure**, same treatment as MusiteDeep's original-repo Python 2.7 blocker.

**Generalized observation, worth recording since it's now happened twice in different forms:** A.1's survey-level check ("weights downloadable: Y") passed both PhosphoSVM and Musite-original, and both failed at actual hands-on use — PhosphoSVM because the download is weights-without-code, Musite-original because the download is a GUI-only installer with no confirmed CLI path (still pending, timeboxed separately). **For older (pre-2015) tools specifically, "weights/model available for download" and "scorable without reimplementation" are not the same claim**, and a survey that only checks for a download link cannot catch this gap — only hands-on invocation does. Worth keeping in mind for Step 4/5 if any further legacy tools are considered.

### DeepPhos (Luo et al.) — **verified runnable, clean pass**

Fetched the repo (`USTC-HIlab/DeepPhos`) and its file tree directly:
- `predict.py` exists, documented as choosing general or kinase-specific prediction via parameters in `main()`.
- Input: CSV (protein name, position, sequence window columns). Output: pipe-separated text (`"Q99440" "3" "0.12992426753"` — protein ID, position, probability).
- **26 real `.h5` weight files confirmed in `models/`** (not empty placeholders — sizes 2.91–3.95MB each), including `model_general_S,T.h5` and `model_general_Y.h5` — exactly the general (non-kinase-specific) models H6 needs.
- Dependencies: `keras==2.0.0`, `tensorflow` backend, `numpy>=1.8.0` — no exact Python/TF version pinned, a 2019-era stack, plausibly needs an older TF/Keras pin to load the `.h5` files correctly (not yet tested end-to-end).
- **No LICENSE file** — run-and-report only, no redistribution, per the existing sign-off (also unrelated to the H14(b) naming-collision requirement, which stays regardless of whether this tool is used).

DeepPhos clears H7 on hands-on inspection, not just the survey. Actual end-to-end run (environment setup, CSV conversion from the frozen sample, inference) is the next concrete step for this tool.

### MusiteDeep — two distinct findings, web server path downgraded, better local path discovered

**(a) The public web server (`musite.net`) is a React single-page app — its actual submission mechanism could NOT be confirmed as scriptable.** The homepage is a JS-rendered SPA (`<div id="root">`, "You need to enable JavaScript to run this app") with no visible form action or API endpoint in the raw HTML. Checked the web server's own source repo (`duolinwang/MusiteDeep_web`) for documented API/curl examples — **none found.** The earlier "10MB/job, 5 concurrent jobs, ~3min/1000 sequences" figures (from the NAR 2020 paper) describe the web UI's *upload limits*, not a confirmed scriptable submission path. **This is now flagged as unconfirmed, not verified viable** — same category of risk PhosphoSVM's CAPTCHA turned out to be, correcting the earlier round's premature "clears the bar" framing.

**(b) Better find: `MusiteDeep_web` repo's own `MusiteDeep/` subfolder is a separate, modern, standalone, scriptable tool — not the same Python 2.7/Theano code as the original 2017 repo.** Fetched its README directly:
- **Python 3.5.2** (tested Ubuntu 16.04, Mac OS High Sierra) — a completely different, far more modern stack than the original repo's Python 2.7/Theano.
- TensorFlow 1.12.0 (or `tensorflow-gpu`), Keras 2.2.4.
- Documented CLI: `python3 predict_multi_batch.py -input [FASTA] -output [prefix] -model-prefix [prefix]`, supports multiple PTM models via semicolon-separated `-model-prefix`.
- **13 pretrained models included directly in `MusiteDeep/models/`** in the repo — no separate weight-hunting needed.
- FASTA in, no CAPTCHA, no browser interaction, no manual step — fully scriptable.

**This was not on the radar in the earlier feasibility round**, which only checked the original `duolinwang/MusiteDeep` repo (Python 2.7) and its `MusiteDeep_Keras2.0` subfolder (confirmed still Python 2.7 underneath) — it never checked the *separate* `MusiteDeep_web` repository's own bundled standalone tool, which turns out to be a different, modern implementation built for the 2020 web-server paper.

**This changes the shape of the MusiteDeep decision** your standing instruction was calibrated to the 2017 Python 2.7/Theano repo's specific difficulty ("do NOT attempt locally... time-box on Kaggle"). This newer repo's constraints (Python 3.5, pip-installable TF 1.12/Keras 2.2.4) are a different and much easier proposition — plausibly runnable via WSL Ubuntu (confirmed present on this machine) without needing Kaggle at all, though not yet attempted or timeboxed. **Flagging this as new information rather than unilaterally acting on it**, since it changes the premise the "don't attempt locally" instruction was based on.

---

**Running count against H7 (need ≥4): PhosphoSVM FAILED, Musite-original PENDING (installer, timeboxed separately), DeepPhos CONFIRMED runnable, PhosF3C checkpoint-verified-live but not yet hands-on tested, MusiteDeep's original web-server path UNCONFIRMED but a better local path (MusiteDeep_web's standalone tool) just surfaced. Not yet at a confirmed 4 — DeepPhos is the only tool that has cleared full hands-on verification so far this round.**

---

## MusiteDeep local path — hands-on verification (2026-09-26)

Per standing instruction, verified before use rather than assumed from the prior round's discovery.

**1. Scriptable entry point — confirmed.** Cloned `duolinwang/MusiteDeep_web` (`study/tools/MusiteDeep_web/`, gitignored). `MusiteDeep/predict_multi_batch.py` (206 lines, read in full) is a real CLI: `argparse`, `-input [FASTA] -output [prefix] -model-prefix [prefix]`, loads Keras weights via `.load_weights()`, writes per-site probability scores to `<output>_results.txt`. Not training code, not a GUI wrapper.

**2. Real weight files — confirmed.** `models/Phosphoserine_Phosphothreonine/` and `models/Phosphotyrosine/` each contain 10-fold CNN ensembles (2,497,500 bytes each) and 10-fold capsule-net ensembles (9,395,096 bytes each) — 40 files total across both PTM types. Verified genuine HDF5 via magic-byte check (`\x89HDF\r\n\x1a\n`) on representative files from both architectures, not empty placeholders or LFS pointers.

**3. End-to-end run — confirmed.** Found a pre-existing WSL Ubuntu conda env (`musitedeep_env`: Python 3.6.13, TensorFlow 1.12.0, Keras 2.2.4 — matches the repo's documented stack). Added missing deps (`scikit-learn==0.24.2`, `pandas==1.1.5`, `pillow`, `h5py==2.10.0` — pinned below 3.x for compatibility with Keras-2.2.4-era `.h5` attribute encoding). Ran `predict_multi_batch.py` against a 1-sequence slice of the repo's own bundled test FASTA (`testdata/Phosphorylation/ST/test_allspecies_sequences.fasta`) with `-model-prefix models/Phosphoserine_Phosphothreonine`: produced varied, non-degenerate per-residue probability scores (e.g. S19: 0.856, T85: 0.072) in `output_results.txt`. Total wall time ~3 minutes, well under the 90-minute timebox that had been set.

**Verdict: MusiteDeep (web-repo `MusiteDeep/` subfolder, WSL `musitedeep_env`) — VERIFIED RUNNABLE.** This supersedes the "web-server path unconfirmed" line above as this tool's execution path — not the original 2017 web-server-first plan, not the Python 2.7/Theano original repo (still infeasible, unchanged). Publication-year note for the leakage-free analysis: this repo's models are a **different, later training run** than the original 2017 tool (the web-repo's own README and model set correspond to the 2020 *NAR* web-server paper) — the correct publication-year cutoff for this specific model version needs confirming against that paper before the leakage-free secondary analysis, not assumed to be 2017.

**Confirmed (2026-09-26):** `MusiteDeep_web`'s own README `### Citation` section names exactly one paper — Wang et al., *Nucleic Acids Research* 48(W1):W140–W146, 2020-07-02 (`doi:10.1093/nar/gkaa275`). **Leakage-free cutoff for this MusiteDeep execution path = 2020**, not 2017.

---

## DeepPhos (Luo et al.) — hands-on end-to-end verification (2026-09-26)

Cloned `USTC-HIlab/DeepPhos` to `study/tools/DeepPhos/` (gitignored). The prior "CONFIRMED runnable" note above was code-inspection only (fetched via GitHub, never actually executed) — running it end-to-end surfaced three genuine Python-2-to-3 porting breaks, all fixed with minimal, mechanical, behavior-preserving edits (no algorithmic or weight changes):

1. `methods/dataprocess_predict.py`: `half_len = (window_size - 1) / 2` produced a `float` under Python 3 (silently used to slice strings, which errors) — changed to `// 2` (integer division; Python 2's `/` on two ints was already floor division, so this restores identical behavior, not a behavior change).
2. `methods/model_n.py`: `from phosnet import Phos` is an implicit-relative import, legal in Python 2 but not Python 3 — changed to the equivalent explicit `from methods.phosnet import Phos`.
3. `methods/phosnet.py`: uses the Keras-1.x functional `merge(list, mode='concat', concat_axis=N)` API at two call sites. Rather than patch call sites (which would require the same downstream numerics to be re-verified), installed the exact pinned `keras==2.0.0` from `requirement.txt` (paired with `tensorflow==1.12.0`, in a new WSL conda env `deepphos_env`, Python 3.6) — Keras 2.0.0 still ships the `keras.legacy` compatibility shims for the old `merge`/`border_mode`/`init`/`bias`/`W_regularizer` argument names (confirmed via the runtime's own deprecation warnings, e.g. *"The `merge` function is deprecated and will be removed after 08/2017... "*), so **no source changes were needed for this one** — installing the tool's own pinned dependency version was sufficient.

**End-to-end run confirmed:** `python predict.py` (default `__main__`, kinase-specific AGC-group model, bundled `test data.csv`) completed without error and reproduced the exact same numbers already recorded in this file's Task C section above (`"Q99440" "3" "0.12992424"` — matches `"0.12992426753"` to displayed precision). This is the first real hands-on execution of DeepPhos this project has done; the general S/T and Y models (`model_general_S,T.h5`, `model_general_Y.h5`) needed for Task B have not yet been run against the frozen sample — that is Task B's own work, not part of this verification.

---

## Task A (final methodology session, 2026-09-26) — PhosF3C and LMPhosSite, timeboxed

### PhosF3C — checkpoint completed, one test prediction VERIFIED, within the 45-minute timebox

Setup began 11:23. `study/tools/phosf3c/PhosF3C-main/model_ckp/` was already partially populated from an earlier session (`conformer_T.pt`, `conformer_Y.pt` complete at 1,212,737,182 bytes each; `lora_S.pt` an incomplete `.part` at 870MB; no `lora_T.pt`/`lora_Y.pt` at all). A full folder-download via `gdown` was started, then redirected to a single targeted file fetch (`gdown <file_id> -O lora_Y.pt`) once the log revealed the true full size (2,615,534,733 bytes — much larger than the conformer files) and that a whole-folder pass would waste time re-touching the two already-complete conformer files. `lora_Y.pt` finished at 11:52 (verified size matches, HDF5-valid).

Ran `predict.py --config_path ./hparams/conformer_Y_cpu.yaml` (device overridden to `cpu`; `type: ['Y']`; checkpoints pointed at the completed Y pair) against the repo's own bundled `data/merged_predict_Y.fasta` (1,200 sequences, 21,798 candidate Y positions) using the existing Windows Python environment (torch 2.10.0+cpu, transformers 5.12.0, peft 0.21.0 — newer than the repo's pinned versions, used as-is since already installed; no compatibility issue observed). ESM2-650M (2.5GB, no cached copy existed anywhere in this project) downloaded automatically via `torch.hub`. Log confirms: LoRA adapter loaded cleanly (`trainable params: 2,729,128 || all params: 653,772,382`), then genuine per-batch prediction progress (28 batches of 64 positions each completed cleanly, ~29s/batch, zero errors) before the run was deliberately stopped — the full 1,200-sequence file was not needed to satisfy "one test prediction works," and letting it finish (~341 batches, ~2.5 hours) was out of scope for a verification step.

**Verdict: PASS, well within the 45-minute timebox** (checkpoint completed ~29 min in, first verified prediction batches by ~37 min in). Runtime extrapolation for Task B: ~29s per 64 positions ⇒ the 10,000-site frozen sample is expected to take roughly 70–80 minutes end-to-end on CPU — plannable, not open-ended. **PhosF3C is locked in for Task B, type Y model only verified this round** (S and T conformer/lora pairs are incomplete — S has a partial lora file, T has no lora file at all; completing those, or restricting Task B's PhosF3C run to Y-type sites only, is a decision for Task B, not assumed here).

### LMPhosSite — NOT COMPLETED, timebox exceeded (process failure, not a tool failure)

Setup began ~11:45 (30-minute timebox, deadline ~12:15). Found the repo already contains real, non-placeholder weights **in-repo** (`models/ST_model.h5`, `models/Y_model.h5`, verified via HDF5 magic-byte check) — better than the prior survey's "Google Drive link, not in-repo" note. `predict.py` is a genuine end-to-end FASTA-in CLI: it runs ProtT5-XL-UniRef50 embedding extraction live via `transformers` (device already hardcoded to `cpu` in the script — no GPU dependency), no separate precompute step required for a basic test. Created a Python 3.8 conda env (`lmphossite_env`, WSL Ubuntu) and installed the pinned core deps (`biopython`, `numpy`, `pandas`, `tensorflow==2.3.1`, `tqdm`, `sentencepiece`) successfully.

**What did not happen:** `torch`/`transformers` were never installed in that env, and `predict.py` was never run. Attention shifted to monitoring PhosF3C's checkpoint download and test run, which overran into LMPhosSite's window; by the time attention returned, wall-clock since LMPhosSite setup began had reached ~49 minutes, already past its 30-minute timebox. Per the standing rule ("do not extend"), stopping here rather than continuing past the box.

**Verdict: NOT RUN this session — timebox exceeded before a test prediction was attempted.** Recorded as a process failure (my attention allocation), not a finding about the tool itself — the bundled real weights and simple CLI make LMPhosSite look more tractable than PhosF3C was, and it remains a documented candidate for a future attempt with its own fresh timebox, not folded into this cycle's tool set.

### Final locked tool set for Task B (3 external tools, short of H7's ≥4 target — recorded as a deviation per instruction)

| Tool | Status | Execution path |
|---|---|---|
| DeepPhos (Luo et al.) | Verified runnable | Local, WSL/Windows |
| MusiteDeep (web-repo `MusiteDeep/` subfolder) | Verified runnable | Local, WSL `musitedeep_env` |
| PhosF3C | Verified runnable (Y-type checkpoint only, this round) | Local CPU, ESM2+LoRA+Conformer |
| LMPhosSite | Not run — timebox exceeded | Not attempted |
| Musite-original, PTransIPs, PhosphoSVM | Not attempted, per instruction | Reasons already on file (installer-only/no CLI, no license/hosting durability, CAPTCHA-gated + weights-without-code respectively) |

**H7 deviation:** 3 of 4 targeted external tools run, not ≥4. G3 remains evaluable (its own criterion requires ≥2 external tools).

---

## Task B — Run and score, final results (2026-09-26, final methodology session)

All three locked tools scored against the frozen 10,000-site sample (checksum re-verified before scoring: `0285523...477a6cb`, match). MusiteDeep and PhosF3C both used a windowed-input optimization instead of full-length sequences — confirmed exact by direct code read (each tool's own window-extraction logic re-derives an identical local window from whatever sequence it is given) and, for PhosF3C, additionally spot-checked empirically (3 sites scored both ways, full-length vs. windowed, matched to ~1e-6). See `build_musitedeep_windows.py` / `build_phosf3c_windows.py` for the full derivation. Full metrics: `study/processed_data/epsd/predictor_audit/task_b_metrics_results.json`. Figure: `task_c_gate_g3_figure.png`.

**Gap direction, stated explicitly and unambiguously (close-out clarification):** `SEI_decile` is assigned by `build_eval_sample.py` via `pd.qcut(SEI_primary, q=10, labels=1..10)`, and `qcut` labels increase with the underlying value — **decile 1 = LOWEST SEI_primary = understudied proteins; decile 10 = HIGHEST SEI_primary = well-studied proteins.** "bottom3" (deciles 1–3) is the understudied end; "top3" (deciles 8–10) is the well-studied end. `gap = AUROC(top3) − AUROC(bottom3)`. **A positive gap means the tool scores better on well-studied proteins than on understudied ones** — the direction the confound thesis predicts.

### Per-decile AUROC, all three tools (direction visible directly, not just as a summary gap)

| Decile (1=understudied → 10=well-studied) | DeepPhos | MusiteDeep | PhosF3C (Y-only, n=1,359 total) |
|---|---|---|---|
| 1 | 0.6902 [0.658,0.723] | 0.6442 [0.612,0.680] | 0.6086 [0.518,0.700] |
| 2 | 0.7226 [0.693,0.751] | 0.6920 [0.658,0.724] | 0.6956 [0.604,0.785] |
| 3 | 0.7450 [0.716,0.776] | 0.7230 [0.692,0.754] | 0.7095 [0.613,0.797] |
| 4 | 0.7437 [0.714,0.773] | 0.7315 [0.701,0.762] | 0.6964 [0.606,0.781] |
| 5 | 0.7685 [0.740,0.796] | 0.7544 [0.724,0.782] | 0.6860 [0.588,0.782] |
| 6 | 0.7537 [0.725,0.784] | 0.7198 [0.689,0.751] | 0.7106 [0.616,0.796] |
| 7 | 0.7638 [0.737,0.792] | 0.7577 [0.728,0.787] | 0.7267 [0.649,0.804] |
| 8 | 0.7738 [0.745,0.803] | 0.7745 [0.748,0.802] | 0.7919 [0.713,0.868] |
| 9 | 0.7408 [0.711,0.772] | 0.7404 [0.711,0.769] | 0.6451 [0.555,0.734] |
| 10 | 0.7527 [0.722,0.781] | 0.7317 [0.702,0.761] | 0.7220 [0.639,0.805] |

All three tools show their lowest AUROC at decile 1 (understudied) and a generally rising trend toward decile 8–10 (well-studied), though non-monotonic in the middle deciles for all three — visible directly in `task_c_gate_g3_figure.png`.

| Tool | n scored / 10,000 | Pooled AUROC | Gap (top3 − bottom3, positive = better on well-studied) |
|---|---|---|---|
| DeepPhos | 10,000 (100%) | 0.7448, CI [0.7358, 0.7541] | **+0.0371**, CI [0.0111, 0.0619] |
| MusiteDeep | 10,000 (100%) | 0.7267, CI [0.7164, 0.7363] | **+0.0632**, CI [0.0350, 0.0876] |
| PhosF3C | 1,359 (13.59%, Y-type checkpoint only — S/T incomplete this round) | 0.6962, CI [0.6685, 0.7255] | **+0.0541**, CI **[-0.0231, 0.1216]** |

**Rank order:** DeepPhos > MusiteDeep > PhosF3C, identical under standard pooled AUROC and SEI-stratified (mean-of-deciles) evaluation. **No change.**

### Leakage-free secondary analysis (Task B.5, pre-registered) — the bottom3-vs-top3 gap, not just pooled AUROC

**Comparability caveat, stated plainly per instruction:** the leakage-free subset restricts POSITIVES to sites first-reported after the model's own publication year; NEGATIVES are the same fixed unlabelled set used in the primary analysis. This means leakage-free pooled AUROC (and gap) is **not directly comparable in absolute terms** to the primary analysis — the positive population differs. The valid comparison is qualitative: does the understudied-vs-well-studied pattern persist, not whether the two numbers sit on the same scale.

| Tool | Leakage-free n positives (of total) | Per-decile n (bottom3 / top3) | Leakage-free pooled AUROC | Leakage-free gap (top3−bottom3) |
|---|---|---|---|---|
| DeepPhos (pub_year=2019) | 1,315 / 5,000 | d1=86, d2=116, d3=138 / d8=146, d9=132, d10=162 (all ≥50, none flagged) | 0.7773, CI [0.7626, 0.7912] | **-0.0031**, CI [-0.0361, 0.0343] |
| MusiteDeep (pub_year=2020) | 1,040 / 5,000 | d1=70, d2=95, d3=114 / d8=123, d9=102, d10=129 (all ≥50, none flagged) | 0.8010, CI [0.7871, 0.8151] | **-0.0173**, CI [-0.0551, 0.0230] |
| PhosF3C (pub_year=2024) | 0 / 515 | n/a | too small to evaluate (n<200) | n/a |

**This is a genuinely important, previously-uncomputed result (Task B.5 was pre-registered but not reported until now): for both DeepPhos and MusiteDeep, the understudied-vs-well-studied gap essentially disappears — and is slightly negative — once restricted to leakage-free positives.** DeepPhos's leakage-free gap (-0.0031, CI includes zero) is a near-total collapse from its already-modest primary gap (+0.0371). MusiteDeep's leakage-free gap (-0.0173, CI includes zero) is not just smaller but reverses sign from its primary gap (+0.0632, the one gap that cleared G3's threshold on its own). **Stated plainly, not smoothed over:** the primary-analysis gap for both tools is at least partly attributable to differential representation of already-well-established (training-era, pre-cutoff) phosphosites across SEI deciles, rather than to a robust generalization gap on genuinely novel, post-training sites specifically. This does not overturn G3's mechanical evaluation (G3's own criterion is defined on the primary/standard analysis, not the leakage-free secondary one — see Task C below), but it materially weakens the strength of the case for *any* of the three tools showing a genuine, training-independent understudied-vs-well-studied gap, and should be reported with equal prominence to the primary-gap numbers in the manuscript, not relegated to a footnote.

## Task C — Gate G3

**Criterion, quoted verbatim (`DECISIONS_HIGH_TIER_PUSH.md` §4, Step 3):**

> **G3 (pre-registered, decides the tier):** Tier A if **≥ 2 external tools** show an understudied-vs-well-studied AUROC gap **≥ 0.05** with 95% CI excluding zero, **or** any change in rank order of tools between standard and SEI-stratified evaluation. Tier B otherwise — the null result is still reported: "existing predictors are robust to study effort" is itself a finding, and it goes in the paper.

**Per-tool evaluation against both conditions of the gap criterion:**

| Tool | Gap ≥ 0.05? | CI excludes zero? | Both conditions met? |
|---|---|---|---|
| DeepPhos | No (0.0371) | Yes | **No** |
| MusiteDeep | **Yes** (0.0632) | **Yes** [0.0350, 0.0876] | **Yes** |
| PhosF3C | Yes (0.0541), but n=1,359/Y-only | **No** — CI [-0.0231, 0.1216] includes zero | **No** |

**Only 1 of 3 tools (MusiteDeep) meets both conditions of the gap criterion — short of the required ≥2.** Rank order is unchanged between standard and SEI-stratified evaluation for all three tools (the second, disjunctive route to Tier A). **Neither route to Tier A is satisfied.**

**G3: TIER B**, evaluated on the pre-registered primary/standard analysis, per the gate's own criterion (the leakage-free analysis is a separate secondary check, not part of G3's decision rule). This is not a clean null on the primary numbers alone: MusiteDeep clears both thresholds on its own, and PhosF3C's point estimate (0.0541) is consistent with a real gap of similar magnitude but lacks the statistical power to confirm it (1,359 Y-type sites, 13.6% of the sample). **However, the leakage-free secondary analysis (see above) cuts the other way**: for the two tools where it could be checked, the gap does not survive restriction to leakage-free positives — it collapses to near-zero (DeepPhos) or reverses (MusiteDeep), both with CIs including zero. **Taken together, the most defensible characterization is: a suggestive primary-analysis gap exists (driven mainly by MusiteDeep), but it does not clearly survive a training-era-independence check, so it should not be reported as a confirmed, robust understudied-vs-well-studied gap.** This is more cautious than either "predictors are robust to study effort" (too strong a null) or "predictors show a confirmed gap" (not supported once leakage-free is considered) — the manuscript should state both results and let the reader see the tension, not resolve it artificially in either direction.

**Caveats attached to this result (per instruction):**

- **Negative-sampling caveat (A.3):** unlabelled positions on low-SEI (understudied) proteins are more likely to contain undiscovered true phosphosites than unlabelled positions on well-studied proteins, which can artificially deflate measured AUROC specifically in the low-SEI deciles. This could inflate the apparent gap seen here (all three tools show lower AUROC in low deciles) or, if the effect is strong enough, partially mask a larger true gap. Direction of bias is not resolvable from this data alone.
- **Training-set recall (corrected — see the leakage-free bottom3-vs-top3 gap section above, not just pooled AUROC):** leakage-free *pooled* AUROC is higher than primary for both DeepPhos (0.745→0.777) and MusiteDeep (0.727→0.801) — both tools generalize fine to demonstrably-post-training sites in the aggregate. **But the understudied-vs-well-studied *gap itself* does not survive leakage-free restriction for either tool** — DeepPhos's collapses to -0.0031 (CI includes zero) and MusiteDeep's reverses to -0.0173 (CI includes zero), down from primary gaps of +0.0371 and +0.0632 respectively. **The correct characterization is that the primary gap is at least partly a training-era-representation artifact, not a robust generalization gap** — this is a materially different and more cautious conclusion than "the gap survives leakage-free" and supersedes that framing. PhosF3C's leakage-free analysis could not run (0 positives post-2024 cutoff, below the pre-registered n≥200 threshold) — expected going in, given the tool's late (Dec. 2024) publication date relative to the corpus.

**Baseline coverage check (performed now — `phaseH1_step3_taskA_final_lock.md` had left this as "not yet performed"):** the legacy H5 embeddings (`legacy/deploy/deepphos-backend/data/sty_embeddings.h5`, 30,355 proteins) overlap with the frozen sample's 8,331 unique proteins at only **2,586/8,331 = 31.04%** — well short of the pre-registered ≥90% bar. **Baselines (PhosConsensus, ESM-2 probe) are skipped, per the pre-registered rule, on a now-verified coverage check rather than an assumption.**

**Draft figure:** `study/processed_data/epsd/predictor_audit/task_c_gate_g3_figure.png` — AUROC by SEI decile, one line per tool, with bootstrap 95% CI error bars. No baselines run this round (see coverage check just above).

---

## PhosF3C S/T completion attempt — pre-registration (2026-09-26, close-out session)

**Committed before any download or run, per explicit instruction:** PhosF3C S/T completion is fulfilment of the already-locked tool set (PhosF3C was locked in Task A; only its Y-type checkpoint pair was complete at that time). This is not a new tool, not a new analysis, and not a re-litigation of G3's criterion. **One timeboxed attempt, 2 hours.** If the S/T checkpoints are not actually published/downloadable within that window, this is recorded as such and the attempt stops — no substitute model, no partial/approximate checkpoint is used in their place. **If the run completes, G3 is re-evaluated exactly once, against its unchanged verbatim criterion, using the fuller PhosF3C result (S/T+Y instead of Y-only), and that outcome is final** — not subject to a further re-run or cherry-pick.

**State before starting, checked directly:**
- `model_ckp/conformer/conformer_T.pt` — present, 1,212,737,182 bytes (byte-identical size to the already-verified `conformer_Y.pt`) — plausibly complete, not yet HDF5/torch-load verified.
- `model_ckp/conformer/conformer_S.pt` — **absent.**
- `model_ckp/lora/lora_T.pt` — **absent.**
- `model_ckp/lora/lora_S.ptx_hy2hl_.part` — **incomplete partial download**, 884,473,856 bytes vs. `lora_Y.pt`'s complete 2,615,534,733 bytes (~34% of expected size, if S is the same size as Y).

Both S and T are currently missing at least one of their two required checkpoint files (a conformer/lora pair is needed together per `hparams/conformer_Y_cpu.yaml`'s structure). Neither can run as-is.

### Outcome: checkpoints downloaded successfully, run exceeded the 2h timebox — stopped, not substituted

**Checkpoints — confirmed available and downloaded.** Listed the Google Drive folder directly (`gdown.download_folder(..., skip_download=True)`): all 6 expected files present (`conformer_{S,T,Y}.pt`, `lora_{S,T,Y}.pt`). Downloaded the 3 missing/incomplete ones — `conformer_S.pt` (1,212,737,182 bytes), `lora_S.pt` and `lora_T.pt` (2,615,534,733 bytes each) — and verified all 6 checkpoint files torch-load cleanly (`torch.load(..., weights_only=True)`, 750/844 tensors as expected). **This part of the pre-registration's contingency ("if not published/downloadable, stop and record") does not apply — they are published and downloadable.**

**Windowed input built** (`build_phosf3c_windows_st.py`, same window-extraction equivalence as the Y-type run): 5,651 S-type + 2,990 T-type sample sites.

**Run launched, both types in parallel (S and T each need a distinct conformer/lora pair, so no cross-contamination risk running concurrently) to make the best use of remaining budget.** Observed batch counts on the actual windowed input: **S = 385 batches** (24,640 candidate positions — S is the most common residue type, so even a ±15-residue window around 5,651 target sites picks up far more incidental S candidates than Y did), **T = 133 batches** (8,512 candidate positions). At the observed first-batch rate (~45–48s/batch, both jobs sharing the same 12 logical cores): **S extrapolates to ≈5.1 hours, T to ≈1.7 hours** — both individually over the pre-registered 2-hour cap, before even accounting for the two jobs competing for CPU.

**Stopped both processes after ~2–3 minutes of actual runtime (1 batch each), per the pre-registered "one timeboxed attempt, 2 hours" rule and the standing "do not extend" instruction.** Real elapsed time from the pre-registration commit (17:34:48) to this stop was ~90 minutes, almost entirely spent on the (necessary, unavoidable) ~4.8GB checkpoint download — leaving no realistic remaining budget for a multi-hour scoring run. **No partial or substitute result is used.** This is a genuine timeout of the *run*, distinct from the *checkpoint-unavailability* contingency the pre-registration anticipated — recorded as such rather than conflated with it.

**Per the pre-registration's own terms, G3 is NOT re-evaluated, since the S/T attempt did not complete.** The Y-only PhosF3C basis (Task B/C above) stands as the final, pre-registered G3 outcome: **Tier B.** PhosF3C's S/T checkpoints are now confirmed available and downloaded for a future session with a realistic multi-hour budget, should the manuscript revision want a fuller PhosF3C result — but that is future work, not part of this cycle's frozen result.
