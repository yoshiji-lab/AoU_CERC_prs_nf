# AoU_CERC_prs_nf — PRS scoring pipeline (Nextflow DSL2)

This repository contains a **local-first Nextflow pipeline** for computing **polygenic risk scores (PRS)** in the **All of Us (AoU) Researcher Workbench** from:

- **PLINK2 pfiles** (`.pgen/.pvar/.psam`) for AoU genotypes, and
- per-trait **GWAS weight files** (usually LD-clumped) containing effect sizes.

The pipeline is optimized for the “**22-job structure**”: **one Nextflow task per chromosome**, where each task loads the chromosome pfile once and scores **all traits** for that chromosome (fast + efficient I/O).

---

## What this pipeline produces

For each trait, you get:

1. **Per-chromosome score outputs** from PLINK2 (`*.sscore`)
2. A **genome-wide PRS** per individual (trait_SUM)
3. A **wide PRS matrix** across all traits (one row per sample, one column per trait)

The pipeline also writes helpful Nextflow diagnostics:

- `pipeline_report.html` (HTML report)
- `pipeline_timeline.html` (Gantt-style timeline)
- `pipeline_trace.txt` (one line per task: status, runtime, etc.)
- `pipeline_dag.svg` (pipeline graph)

---

## Repository layout

```
AoU_CERC_prs_nf/
├── main.nf
├── nextflow.config
├── nextflow                # local Nextflow launcher (optional)
├── modules/
│   ├── score_trait_chr.nf  # legacy: 1 task per (trait × chr)
│   ├── merge_trait.nf      # merge chr-level sscore → genome-wide per trait
│   └── build_matrix.nf     # merge all traits → wide matrix
├── bin/
│   ├── make_scorefile.py       # weights + PVAR → PLINK2 score file
│   ├── merge_sscore_parts.py   # sum chr scores → genome-wide PRS per trait
│   └── build_prs_matrix.py     # merge per-trait PRS → wide matrix
└── docker/
    └── Dockerfile
```

### Key files (what they do)

**`main.nf`**
- Orchestrates the full workflow:
  - discovers weight files (one per trait)
  - iterates over chromosomes
  - runs scoring, merging, and matrix building

**`nextflow.config`**
- Local-only defaults for AoU Workbench (safe for shared GCS buckets)
- Sets:
  - local `outdir` and `workDir`
  - `stageInMode = 'symlink'` (avoids copying huge `.pgen` files)
  - trace/report/timeline/dag outputs

**`modules/score_trait_chr.nf` (legacy mode)**
- Runs scoring for **one trait on one chromosome**.
- Internally:
  - builds a PLINK2 score file via `bin/make_scorefile.py`
  - runs `plink2 --score …` to produce `*.sscore`

> **Note:** The refactored “one task per chromosome scoring all traits” module is typically named `modules/score_chr_all_traits.nf` and is referenced by `main.nf` in the 22-job design.
> If your checkout doesn’t contain it yet, you can either add it or temporarily switch `main.nf` back to the legacy `SCORE_TRAIT_CHR` approach.

**`modules/merge_trait.nf`**
- Groups all `*.sscore` parts for a trait and merges them into **one genome-wide PRS file**:
  - `TRAIT.genomewide.sscore.tsv`

**`modules/build_matrix.nf`**
- Merges all per-trait genome-wide files into:
  - `prs_matrix_ALL_TRAITS.tsv`

**`bin/make_scorefile.py`**
- Converts GWAS weights + PVAR into a PLINK2 score file:
  - merges by `(chromosome, base_pair_location)`
  - handles direct vs swapped allele alignment
  - optional strand-flip (`--allow_strand_flip`)

**`bin/merge_sscore_parts.py`**
- Reads all chromosome `*.sscore` files for a trait
- Sums per-chromosome scores into a single `TRAIT_SUM`

**`bin/build_prs_matrix.py`**
- Produces one wide table with `FID IID` + one PRS column per trait

---

## Requirements

### Tools
- Nextflow (recommended: **24.10.4** as used in AoU training)
- PLINK2 available on PATH as `plink2`
- Python 3 with:
  - `pandas`
  - `numpy`

### Inputs
You need two sets of local files:

1) **Genotype pfiles** (one per chromosome)
```
<pfile_dir>/snpqc.chr1.pgen
<pfile_dir>/snpqc.chr1.pvar
<pfile_dir>/snpqc.chr1.psam
...
<pfile_dir>/snpqc.chr22.*
```

2) **Weights files** (one file per trait)
- Globbed by `--weights_glob` (e.g. `.../weights/*.tsv.gz`)
- Must include columns:
  - `chromosome`
  - `base_pair_location`
  - `effect_allele`
  - `other_allele`
  - `beta`

Trait names are derived from filenames, with these suffixes stripped:
- `.tsv.gz`, `.csv.gz`, `.tsv`, `.csv`, `.gz`
- and also a trailing `_clumped`

---

## Update Nextflow (AoU training commands)

AoU Workbench often has a preinstalled Nextflow. First check what you have:

```bash
nextflow -v
```

Install the recommended version (24.10.4) in the current directory:

```bash
curl -s -L https://github.com/nextflow-io/nextflow/releases/download/v24.10.4/nextflow | bash

# confirm
./nextflow -v
```

✅ **From here on, always run with the local Nextflow binary**:

```bash
./nextflow run main.nf ...
```

---

## Quickstart (local AoU Workbench run)

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

**Run only selected traits**:

```bash
./nextflow run main.nf -c nextflow.config \
  --traits "VATadj,BMIadj" \
  --chroms 21..22 \
  -resume
```

**Full autosomes**:

```bash
./nextflow run main.nf -c nextflow.config --chroms 1..22 -resume
```

**Increase parallelism on a bigger VM** (careful with RAM/CPU):

```bash
./nextflow run main.nf -c nextflow.config --max_forks 4 -resume
```

**Enable strand-flip matching** (only if you truly need it):

```bash
./nextflow run main.nf -c nextflow.config --allow_strand_flip true -resume
```

---

## Run inside a tmux session (recommended)

Using tmux ensures your run **survives browser disconnects, VPN drops, or idle timeouts**.

### 1) Start a new tmux session

```bash
tmux new -s prs_nf
```

### 2) Run the pipeline *inside tmux*

```bash
cd /home/jupyter/workspaces/allofusgwasonmetabolictraits/prs_nf/AoU_CERC_prs_nf

./nextflow run main.nf \
  -c nextflow.config \
  --chroms 21..22 \
  -resume \
  2>&1 | tee run.nextflow.stdout.log
```

### 3) Detach (leave it running)
Press:

- **Ctrl + b**, then **d**

### 4) Reattach later

```bash
tmux attach -t prs_nf
```

### 5) See running sessions / kill session

```bash
tmux ls

# kill the session when done
 tmux kill-session -t prs_nf
```

✅ Tip: Your most important live debug file is:

```bash
tail -f .nextflow.log
```

---

## Outputs

All outputs are written under `--outdir`.

### PRS results

- `outdir/per_chr/`
  - per-trait × per-chr scoring outputs (`*.sscore`) and logs

- `outdir/prs_merged/`
  - `TRAIT.genomewide.sscore.tsv` (FID, IID, TRAIT_SUM)
  - `prs_matrix_ALL_TRAITS.tsv` (FID, IID + one column per trait)

### Nextflow diagnostics

- `outdir/pipeline_trace.txt`
- `outdir/pipeline_report.html`
- `outdir/pipeline_timeline.html`
- `outdir/pipeline_dag.svg`

---

## How to read the Nextflow HTML report & logs

### Viewing the HTML outputs in AoU Workbench

The HTML files are written under your `--outdir` (e.g. `pipeline_report.html`). In the AoU Workbench Jupyter UI you can:

1. Open the File Browser (left sidebar)
2. Navigate to your `outdir` folder
3. Click `pipeline_report.html` or `pipeline_timeline.html` to open it in a new tab

If your browser does not render it directly, you can download the file from the Jupyter UI and open it locally in Chrome/Firefox.


### 1) `pipeline_report.html`
This is your “what happened?” report.

What to look for:
- **Processes tab**: which step ran (scoring / merging / matrix)
- **Task durations**: slow chromosomes/traits
- **Failures**: which task failed and exit code

### 2) `pipeline_timeline.html`
This is your “parallelization” report.

Look for:
- whether chromosome tasks ran concurrently as expected
- long-running stragglers

### 3) `pipeline_trace.txt`
This is a machine-readable log per task.

Helpful columns:
- `process`, `tag`, `status`, `exit`
- `start`, `complete`, `duration`
- `cpus`, `memory`

### 4) `.nextflow.log`
This is the global Nextflow runtime log.

First place to check when something crashes:

```bash
less .nextflow.log
```

### 5) Per-task work directories
Nextflow stores each task execution under `work/`.

For any failing task, you can inspect:
- `.command.sh` (exact executed script)
- `.command.out` and `.command.err`

Example:

```bash
ls -lh work/
cat work/<hash>/.command.err
```

---

## How `-resume` works (and how to not lose progress)

Nextflow caches tasks based on:
- the process code
- the exact inputs
- the parameters

When you rerun **with the same `workDir`** and add `-resume`, Nextflow will:
- **skip** tasks that completed successfully previously
- **re-run only** tasks that failed or changed

✅ Best practice:
- Keep `--workdir` stable for a project run
- Re-run with the same command + `-resume`

Example after a VM disconnect:

```bash
./nextflow run main.nf -c nextflow.config --chroms 1..22 -resume
```

---

## Troubleshooting

### “It says it can’t find plink2”
Make sure PLINK2 is on PATH:

```bash
which plink2
plink2 --version
```

Or point to it explicitly:

```bash
./nextflow run main.nf -c nextflow.config --plink2_bin /path/to/plink2 -resume
```

### “I got an empty score file / no overlaps”
This usually means:
- your weights don’t overlap the genotype build, or
- the trait was on a different reference/build, or
- alleles didn’t match and `--allow_strand_flip` is off

Look at:
- `outdir/per_chr/*.scoring.log`

### “I changed something and it recomputed everything”
That’s expected if:
- you changed the pipeline code, or
- you changed parameters / inputs,

because task hashes change and Nextflow treats tasks as new.

---

## (Optional) Containerized execution

A minimal container is provided in `docker/Dockerfile` with PLINK2 + pandas/numpy.

If you want Nextflow to run processes inside the container, you can extend the config with:

```groovy
process.container = '<your-built-image>'
```

(Only do this if you need strict reproducibility across environments.)

---

## Citation / acknowledgement

This pipeline wraps PLINK2 scoring and reproduces the PRS aggregation logic used in the companion notebook (allele-aware scorefile creation → per-chr scoring → genome-wide merging → wide matrix).
