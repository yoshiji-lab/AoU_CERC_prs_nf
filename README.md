# AoU_CERC_prs_nf — PRS scoring pipeline (Nextflow DSL2, AoU Workbench)

This repository contains a **local-first Nextflow pipeline** for computing **polygenic risk scores (PRS)** in the **All of Us (AoU) Researcher Workbench**, using:

- **AoU PLINK2 pfiles** (`.pgen/.pvar/.psam`) per chromosome, and  
- per-trait **GWAS weight files** (typically LD-clumped) with effect sizes

✅ **Key design choice (“22-job structure”)**  
This pipeline runs **one job per chromosome**, and inside that chromosome job it scores **all traits**.  
That means each chromosome’s `.pgen/.pvar/.psam` is “loaded” once and reused, which is **much faster and more I/O efficient** than running `(trait × chr)` as separate jobs.

---

## What this pipeline produces

### For each trait
1. **Per-chromosome PLINK2 scoring outputs** (`*.sscore`)
2. A **genome-wide PRS per individual** (a single “SUM” score per sample)

### Across all traits
3. A **wide PRS matrix** (one row per sample, one column per trait)

### Nextflow diagnostics
The pipeline writes standard Nextflow run artifacts into your `--outdir`:

- `pipeline_report.html` (HTML report)
- `pipeline_timeline.html` (timeline / concurrency view)
- `pipeline_trace.txt` (task table)
- `pipeline_dag.svg` (pipeline DAG)

---

## Repository layout (core files)

```
.
├── main.nf
├── nextflow.config
├── nextflow-1.config          # optional alt config (older / different inputs)
├── nextflow                  # local Nextflow launcher (recommended to use)
└── README_run.md             # legacy notes (Batch/GCS style; not recommended)
```

### What each file does

#### `main.nf`
The main Nextflow workflow:
- parses chromosomes (`--chroms`)
- discovers weight files (`--weights_glob`)
- optionally filters traits (`--traits`)
- builds a single `weights.list` file (absolute file paths)
- runs **one job per chromosome** via `SCORE_CHR_ALL_TRAITS`
- groups `*.sscore` by trait and merges them (`MERGE_TRAIT`)
- builds final wide matrix (`BUILD_MATRIX`)

#### `nextflow.config`
Your **local-only safe config** for AoU Workbench.
- avoids Google Batch + avoids `gs://` working dirs by default
- uses `stageInMode = 'symlink'` (important: avoids copying giant pfiles)
- emits trace/report/timeline/dag into `--outdir`
- controls concurrency with `--max_forks`

#### `nextflow` (local binary)
This is the **recommended Nextflow executable** (pinned to the version you installed).
Run the pipeline with:

```bash
./nextflow run main.nf ...
```

#### `README_run.md`
Older “Batch/GCS-style” examples.  
**Not recommended** if you want to avoid shared-bucket overwrite risk.

---

## Requirements

### Tools
- **Nextflow** (recommended: **24.10.4**, from AoU training)
- **PLINK2** available as `plink2` on PATH
- Local filesystem space for `workdir` and outputs

### Inputs

#### 1) Genotype pfiles (one per chromosome)
Example layout:

```
<pfile_dir>/snpqc.chr1.pgen
<pfile_dir>/snpqc.chr1.pvar
<pfile_dir>/snpqc.chr1.psam
...
<pfile_dir>/snpqc.chr22.*
```

You control naming via:

- `--pfile_dir`
- `--pfile_pref` (must contain `{chrom}`, e.g. `snpqc.chr{chrom}`)

#### 2) Weight files (one per trait)
Weight files are discovered by:

- `--weights_glob "/path/to/weights/*.tsv.gz"`

Trait names are derived from filenames by stripping:
- `.tsv.gz`, `.csv.gz`, `.tsv`, `.csv`, `.gz`
- and a trailing `_clumped`

**Expected columns in weights** (minimum):
- `chromosome`
- `base_pair_location`
- `effect_allele`
- `other_allele`
- `beta`

---

## Update Nextflow (AoU training commands)

AoU Workbench sometimes has a preinstalled Nextflow. Check it:

```bash
nextflow -v
```

Install **Nextflow 24.10.4** in the current directory:

```bash
curl -s -L https://github.com/nextflow-io/nextflow/releases/download/v24.10.4/nextflow | bash
chmod +x nextflow
./nextflow -v
```

✅ From here on, **always run with the local binary**:

```bash
./nextflow run main.nf ...
```

---

## Quickstart (recommended local AoU Workbench run)

From the repo root:

```bash
./nextflow run main.nf \
  -c nextflow.config \
  --chroms 21..22 \
  --weights_glob "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/inputs/test/*.tsv.gz" \
  --pfile_dir "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/pgen_dir" \
  --pfile_pref "snpqc.chr{chrom}" \
  --outdir "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/prs_out_test" \
  --workdir "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/nf_work/prs_out_test" \
  -resume
```

### Useful options

Run only selected traits:

```bash
./nextflow run main.nf -c nextflow.config \
  --traits "VATadj,BMIadj" \
  --chroms 21..22 \
  -resume
```

All autosomes:

```bash
./nextflow run main.nf -c nextflow.config --chroms 1..22 -resume
```

Control max parallel chromosome jobs (local VM safety):

```bash
./nextflow run main.nf -c nextflow.config --max_forks 2 -resume
```

Enable strand-flip matching (only if needed):

```bash
./nextflow run main.nf -c nextflow.config --allow_strand_flip true -resume
```

---

# Running safely inside tmux (so it never dies)

If your browser disconnects, VPN drops, or the AoU session times out, `tmux` keeps the run alive.

### 1) Start (or reuse) a tmux session

```bash
tmux new -A -s prs_nf
```

### 2) Run the pipeline *inside tmux* and log stdout

```bash
cd /home/jupyter/workspaces/allofusgwasonmetabolictraits/prs_nf/AoU_CERC_prs_nf

./nextflow run main.nf -c nextflow.config --chroms 21..22 -resume \
  2>&1 | tee run.nextflow.stdout.log
```

### 3) Detach and leave it running

Press:
- `Ctrl + b` then `d`

### 4) Reattach later

```bash
tmux attach -t prs_nf
```

### 5) Check sessions / kill when done

```bash
tmux ls
tmux kill-session -t prs_nf
```

### Live debugging while it runs

The **best live Nextflow log** is always:

```bash
tail -f .nextflow.log
```

---

## Output structure

All pipeline outputs go under `--outdir`.

Typical structure:

- `outdir/per_chr/`
  - per-chromosome scoring outputs (each trait produces `*.sscore`)

- `outdir/prs_merged/`
  - per-trait genome-wide PRS files
  - final matrix across all traits

- `outdir/pipeline_trace.txt`
- `outdir/pipeline_report.html`
- `outdir/pipeline_timeline.html`
- `outdir/pipeline_dag.svg`

---

# How to read the HTML report & logs (what to click / what it means)

## Viewing HTML in AoU Workbench
In the AoU Jupyter file browser:
1. navigate to your `--outdir`
2. click `pipeline_report.html` or `pipeline_timeline.html`
3. open in a new browser tab

If it doesn’t render, download it and open locally in Chrome.

---

## 1) `pipeline_report.html` (the “what happened?” report)
This is your best “one screen summary”.

Look for:
- **Processes section**: which processes ran (scoring / merge / matrix)
- **Task table**: runtimes and failures
- **Failures tab**: which task failed + exit code

Common use:
- “Which chromosome was slow?”
- “Which step failed first?”

---

## 2) `pipeline_timeline.html` (the “parallelization” report)
This shows **concurrency over time** (Gantt chart style).

Look for:
- Are chromosome tasks running in parallel?
- Is one chromosome a straggler (usually chr1/chr2)?
- Are you under-utilizing CPU (max_forks too low)?

---

## 3) `pipeline_trace.txt` (the “spreadsheet” log)
A tabular run record per task.

Useful columns:
- `process`, `tag`, `status`, `exit`
- `start`, `complete`, `duration`
- `cpus`, `memory`
- `hash` (unique task ID)

This is great for quick grep/filters:

```bash
column -t outdir/pipeline_trace.txt | less -S
```

---

## 4) `.nextflow.log` (the global run log)
This is the **first file to check if anything crashes**.

```bash
less .nextflow.log
```

What it tells you:
- why Nextflow stopped
- which process crashed
- stack traces for session/lock problems

---

## 5) Work directory debugging (deep dive)
Every task runs in a unique directory under `work/`.

Inside a failing task dir you’ll see:
- `.command.sh`  (exact script that ran)
- `.command.out` (stdout)
- `.command.err` (stderr)
- `.exitcode`

How to inspect:
```bash
ls -lh work/
cat work/<hash>/.command.err
cat work/<hash>/.command.sh
```

This is the **most foolproof** way to debug *exactly what ran*.

---

# How `-resume` works (and how to not lose progress)

Nextflow **caches each task** based on:
- the process code
- the input files
- the parameters

When you re-run with the **same `--workdir`** and add `-resume`:

✅ tasks that already completed successfully are **skipped**  
✅ only missing/failed/changed tasks are **recomputed**

### Best practices for resume
- keep `--workdir` stable for the whole project run
- always rerun with the same command + `-resume`

Example:
```bash
./nextflow run main.nf -c nextflow.config --chroms 1..22 -resume
```

### When will it recompute things?
It’s normal for Nextflow to rerun tasks if you changed:
- pipeline code (main.nf / modules)
- parameter values
- input weights
- genotype files

Because those changes alter the task hash.

---

# Troubleshooting

## “Unable to acquire lock on session …”
This happens if:
- you actually have a pipeline still running, OR
- a prior run was killed abruptly and left a stale lock

First check if another run is active:
```bash
ps aux | grep nextflow
```

If nothing is running and you’re sure it’s stale, the fix is usually:
- rerun with a *new* `-name`, OR
- remove the stale session cache directory mentioned in the error

(Always prefer the least destructive fix: new run name.)

---

## “plink2 not found”
Check:
```bash
which plink2
plink2 --version
```

Or pass a custom binary path:
```bash
./nextflow run main.nf -c nextflow.config --plink2_bin "/path/to/plink2" -resume
```

---

## “No overlap / empty results”
Usually:
- weights are different build than genotypes
- alleles don’t match and strand-flip is off
- chrom/pos IDs don’t match

What to inspect:
- per-chrom scoring logs under `outdir/per_chr/`
- `.command.err` for the failing task

---

# Presentation / Talk Track (speaker notes)

Use this section as a **script** when presenting the pipeline.

## 30-second overview (what the pipeline is)
> “This Nextflow DSL2 pipeline computes polygenic risk scores in All of Us using PLINK2 pfiles and per-trait GWAS weight files.  
> The key optimization is we run **one job per chromosome**, and score **all traits inside that job**, so we don’t repeatedly reload the same large genotype files.  
> Outputs include per-chromosome `.sscore` files, genome-wide PRS per trait, and a final PRS matrix across all traits.”

---

## What problem it solves (why we built it)
- We want PRS for many traits (dozens to hundreds)
- Genotypes are huge → I/O dominates runtime if you score `(trait × chr)`
- AoU Workbench environment can disconnect → need robustness
- Need reproducible outputs + resumability

---

## Key design decisions (what makes it efficient)
1. **One chromosome task scores all traits**
   - minimizes genotype file staging / repeated reads
2. **Local-only config by default**
   - avoids shared `gs://` overwrite risk in organization buckets
3. **Symlink staging for pfiles**
   - avoids copying massive `.pgen` into every task directory
4. **Nextflow caching (`-resume`)**
   - restartable after disconnects or failures

---

## Walkthrough: what each step does (pipeline flow)
1. **Discover weights**
   - glob all weights files from `--weights_glob`
   - derive trait names from filenames
2. **Build chromosome job inputs**
   - for each chromosome → load matching `.pgen/.pvar/.psam`
3. **Scoring**
   - for each chromosome job → run PLINK2 scoring for all traits
   - produces `trait.chrN.scored.sscore`
4. **Merge per trait**
   - group `*.sscore` across chromosomes for a trait
   - sum to genome-wide PRS per individual
5. **Build matrix**
   - combine all trait genome-wide scores into one matrix

---

## How to “prove it ran correctly” (demo checklist)
- Check `pipeline_report.html` → all tasks succeeded
- Check `pipeline_trace.txt` → no failed tasks (`exit != 0`)
- Confirm outputs exist:
  - `outdir/per_chr/*.sscore`
  - `outdir/prs_merged/` final trait PRS
  - final PRS matrix across traits

---

## How to read the HTML outputs (what to say while showing them)
### `pipeline_report.html`
- “This is the best summary page: processes + runtimes + failures.”
- “If something fails, I immediately use this to identify the exact task.”

### `pipeline_timeline.html`
- “This shows parallelization — I expect multiple chromosomes running at once.”
- “Long bars highlight bottlenecks like chr1/chr2.”

### `pipeline_trace.txt`
- “This is the machine-readable task table — great for quick performance debugging.”

---

## Debugging story (what you check first when it fails)
1. `.nextflow.log` for the overall crash reason
2. failing task work directory:
   - `.command.err` for the real error message
   - `.command.sh` to rerun manually if needed

---

## Resume story (why -resume matters)
> “If my VM disconnects, I rerun the same command with `-resume`.  
> Nextflow skips anything that already finished successfully and only recomputes missing pieces.  
> As long as `--workdir` is unchanged, we don’t lose progress.”

---

## Closing (what this enables downstream)
- PRS matrix can feed:
  - association models
  - clustering/phenotyping
  - correlation with IDPs / metabolic traits
  - stratified PRS analyses by ancestry/sex

---

## Recommended “demo run” commands for a live presentation

Small test run:
```bash
./nextflow run main.nf -c nextflow.config --chroms 21..22 --traits "VATadj,BMIadj" -resume
```

Full autosomes:
```bash
./nextflow run main.nf -c nextflow.config --chroms 1..22 -resume
```

Inside tmux:
```bash
tmux new -A -s prs_nf
./nextflow run main.nf -c nextflow.config --chroms 1..22 -resume 2>&1 | tee run.stdout.log
```

