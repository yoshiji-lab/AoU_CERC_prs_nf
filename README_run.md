# PRS Nextflow pipeline (AoU / Google Batch)

## Run
nextflow run main.nf \
  -c nextflow.config \
  --weights_glob "gs://.../*.tsv.gz" \
  --pfile_dir "gs://.../pgen_dir" \
  --pfile_pref "snpqc.chr{chrom}" \
  --outdir "gs://.../prs_out" \
  --workdir "gs://.../nf_work/prs" \
  --threads 8 \
  --memory "16 GB" \
  --disk "200 GB" \
  --max_forks 50 \
  -resume

## Outputs
- gs://.../prs_out/per_chr/     (per trait×chr score + sscore)
- gs://.../prs_out/prs_merged/  (per-trait genomewide + wide matrix)
- gs://.../prs_out/trace.txt, report.html, timeline.html, dag.svg

## Resume after VM death
As long as --workdir points to the same gs:// location, rerun with -resume.

