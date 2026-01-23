# AoU_CERC_prs_nf (Nextflow) — PRS scoring in **22 chromosome jobs** (recommended)

This repository computes **polygenic risk scores (PRS)** in the **All of Us (AoU) Researcher Workbench** using:
- AoU genotype files in **PLINK2 pfile** format (`.pgen/.pvar/.psam`), one set per chromosome, and
- GWAS **weight files** (usually LD-clumped), one file per trait.

✅ **This README explains only the recommended way to run the pipeline:**  
**the “22-job structure” = 1 job per chromosome (chr1…chr22)**.

---

## The one idea you need to understand: “22-job structure”

If you are new to PRS scoring, the slow part is often not the scoring itself—it's repeatedly opening huge genotype files and doing setup work (like allele frequency calculations).

### What “22-job structure” means (plain English)

- Human autosomal DNA is split into **22 chromosomes** (chr1–chr22).
- AoU genotype data is stored as **one big file-set per chromosome**:
  - `chr7.pgen`, `chr7.pvar`, `chr7.psam`  (and similarly for each chromosome)
- If you score **one trait at a time**, you end up doing this many times:
  - “Open chr7 files → do setup → score trait A”  
  - “Open chr7 files → do setup → score trait B”  
  - “Open chr7 files → do setup → score trait C”  
  …which is slow and costs more compute.

### What this pipeline does instead

For each chromosome **once**:

1. **Load the chromosome pfiles once** (e.g., chr7)
2. **Score ALL traits for that chromosome in one PLINK2 run**
3. Save per-trait per-chromosome results
4. After all chromosomes finish, merge them to build final PRSs

So the pipeline runs **22 chromosome jobs** (chr1…chr22).  
Each chromosome job runs **one PLINK2 scoring pass** for all traits.

That’s the core optimization.

---

## What you will get at the end

You will get:

1. **Per-chromosome PRS outputs** (for debugging + merging)
2. **Per-trait genome-wide PRS** (one PRS per person per trait)
3. **A wide PRS matrix** (one row per person, one PRS column per trait)

And Nextflow run reports (HTML) so you can see what ran and how long it took.

---

## Folder structure (what matters)

```
AoU_CERC_prs_nf/
├── main.nf
├── nextflow.config
├── nextflow                      # local Nextflow launcher (recommended)
├── modules/
│   ├── score_chr_all_traits.nf    # ✅ core: 1 job per chromosome, score ALL traits
│   ├── merge_trait.nf             # merge chr-level scores → genome-wide PRS per trait
│   ├── build_matrix.nf            # merge all traits → final wide matrix
│   └── score_trait_chr.nf         # 🚫 IGNORE for now (older structure, not recommended)
├── bin/
│   ├── make_scorefile.py          # builds PLINK2 score file for a trait+chrom
│   ├── split_multisscore.py       # splits combined .sscore → per-trait .sscore
│   ├── merge_sscore_parts.py      # merges chr scores → genome-wide trait PRS
│   └── build_prs_matrix.py        # builds the final wide PRS matrix
└── docker/
    └── Dockerfile
```

**Important:** `modules/score_trait_chr.nf` exists in the repo but **is not recommended** and **is not explained here**.

---

## Before you run: what you need

### 1) Nextflow
AoU often has a preinstalled Nextflow, but we recommend using a **local pinned version** (same as AoU training).

Check version:
```bash
nextflow -v
```

Install Nextflow `v24.10.4` locally inside this repo folder:
```bash
curl -s -L https://github.com/nextflow-io/nextflow/releases/download/v24.10.4/nextflow | bash
chmod +x ./nextflow
./nextflow -v
```

From now on, run Nextflow as:
```bash
./nextflow run main.nf ...
```

### 2) PLINK2
You need `plink2` available (or provide its path via `--plink2_bin`).

Check:
```bash
plink2 --version
```

### 3) Inputs
You need two inputs:

#### (A) Genotypes: pfiles per chromosome
You need files like:
- `snpqc.chr1.pgen`, `snpqc.chr1.pvar`, `snpqc.chr1.psam`
- ...
- `snpqc.chr22.pgen`, `snpqc.chr22.pvar`, `snpqc.chr22.psam`

You will tell the pipeline where they are using:
- `--pfile_dir` (folder)
- `--pfile_pref` (prefix template, e.g. `snpqc.chr{chrom}`)

#### (B) GWAS weights (one file per trait)
These are your clumped GWAS files. Example filenames:
- `VAT_to_Android_clumped.tsv.gz`
- `Leg_lean_clumped.tsv.gz`
- etc.

You will tell the pipeline which files to use with:
- `--weights_glob "/path/to/weights/*.tsv.gz"`
- optionally `--traits "Trait1,Trait2"` to run a subset

Your weights must contain columns like:
- `chromosome`
- `base_pair_location`
- `effect_allele`
- `other_allele`
- `beta`

---

## Recommended way to run (always): use tmux

AoU notebooks/SSH sessions can disconnect. `tmux` keeps your run alive.

Start a tmux session:
```bash
tmux new -s prs_nf
```

Run your pipeline inside it:
```bash
cd /home/jupyter/workspaces/allofusgwasonmetabolictraits/prs_nf/AoU_CERC_prs_nf
./nextflow run main.nf -c nextflow.config -resume 2>&1 | tee run.stdout.log
```

Detach (leave it running):
- Press `Ctrl+b`, then `d`

Reattach later:
```bash
tmux attach -t prs_nf
```

---

## Quick test run (recommended first time)

Run only **chr21–22** and only **2 traits** to verify everything works.

Example:
```bash
./nextflow run main.nf -c nextflow.config -resume   --chroms 21..22   --traits "VAT_to_Android,Leg_lean"   --weights_glob "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/inputs/00.Adiposity/*.tsv.gz"   --pfile_dir  "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/pgen_dir"   --pfile_pref "snpqc.chr{chrom}"   --outdir  "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/prs_out_test"   --workdir "/home/jupyter/workspaces/allofusgwasonmetabolictraits/urvashi_analysis/nf_work/prs_out_test"
```

What success looks like:
- you see **one job per chromosome**
- inside each chromosome job, you see **one PLINK2 scoring run** (not one per trait)

---

## Full run (chr1–22, all traits)

```bash
./nextflow run main.nf -c nextflow.config -resume   --chroms 1..22
```

---

## Key parameters (simple explanations)

- `--chroms`  
  Which chromosomes to run. Examples:
  - `21..22`
  - `1..22`
  - `1,2,7`

- `--weights_glob`  
  A file pattern that matches your weight files.
  Example: `"/path/to/weights/*.tsv.gz"`

- `--traits` *(optional)*  
  Comma-separated list of trait names to run (for testing).
  Example: `"VAT_to_Android,Leg_lean"`  
  If omitted, pipeline runs **all weights matched by `--weights_glob`**.

- `--pfile_dir`  
  Folder containing chromosome pfiles.

- `--pfile_pref`  
  Prefix template for each chromosome. Must include `{chrom}`.
  Example: `snpqc.chr{chrom}`

- `--outdir`  
  Where final outputs + reports go.

- `--workdir`  
  Nextflow cache/work directory. Keep it stable so `-resume` works well.

- `--plink2_bin` *(optional)*  
  Path to plink2 if not on PATH.

- `--allow_strand_flip` *(optional, default false)*  
  Only enable if you are sure you need it (can hide mismatches).

---

## Where outputs go

Inside `--outdir` you will typically see:

### A) Per-chromosome outputs (intermediate)
- `per_chr/TRAIT.chrX.scored.sscore`

These are useful for debugging *which variants scored* and *how many overlapped*.

### B) Genome-wide PRS per trait
- `prs_merged/TRAIT.genomewide.sscore.tsv`

This is one PRS column per trait per person.

### C) Final PRS matrix (all traits together)
- `prs_merged/prs_matrix_ALL_TRAITS.tsv`

One row per person, one PRS column per trait.

### D) Nextflow reports (great for first-time users)
- `pipeline_report.html` (summary of tasks, failures, durations)
- `pipeline_timeline.html` (timeline view)
- `pipeline_trace.txt` (tabular trace)
- `pipeline_dag.svg` (pipeline DAG)

Open the HTML files in the AoU notebook file browser or download them locally.

---

## How to tell the pipeline is doing the optimized “load once per chromosome” scoring

For each chromosome, look at the chromosome scoring log (published under `per_chr/` or process logs):

You want to see:
- multiple traits being prepared (`make_scorefile` messages)
- then **a single PLINK2 scoring run** for that chromosome

If you see “Calculating allele frequencies…” repeating 30+ times for the same chromosome, that means PLINK2 is being launched per trait (not optimized). This pipeline is designed to avoid that.

---

## If something fails: how to debug

### 1) Look at `.nextflow.log`
```bash
tail -n 200 .nextflow.log
```

### 2) Find the work directory for the failed task
Nextflow prints something like:
```
Work dir: /.../nf_work/.../ab/cdef123...
```
Go there and inspect:
- `.command.sh` (exact command that ran)
- `.command.err` (errors)
- `.command.out` (stdout)

### 3) Rerun with `-resume`
Once fixed:
```bash
./nextflow run main.nf -c nextflow.config -resume ...
```

---

## Common problems (and fixes)

### “plink2: command not found”
Either add plink2 to PATH or pass:
```bash
--plink2_bin /path/to/plink2
```

### “Empty score file / no overlaps”
That trait had no variants overlapping that chromosome (or alleles didn’t match).
This is not always an error—some traits have very few lead SNPs on some chromosomes.

### “unrecognized arguments: --score_name”
Your `bin/make_scorefile.py` must support the `--score_name` argument.  
Update `make_scorefile.py` to include it (this pipeline relies on it to name score columns correctly).

---

## Notes for collaborators (what you should *not* do)

- Do **not** run `score_trait_chr.nf` unless you are explicitly asked to.
- Do **not** run “one trait at a time” scoring loops — it is slower and more expensive.
- Do **not** override configs to use shared org-wide GCS paths unless you know what you're doing.
  This repo is intended to be run **locally in your workspace** to avoid accidental overwrites.

---

## Contact / questions

If you are running this for the first time and something is confusing, the fastest way to debug is to share:
- the command you ran,
- the last ~50 lines of `.nextflow.log`, and
- the failed task’s `.command.err`.

