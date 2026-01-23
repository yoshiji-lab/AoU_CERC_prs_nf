nextflow.enable.dsl = 2

process SCORE_CHR_ALL_TRAITS {
  tag { "chr${chrom}" }

  publishDir "${params.outdir}/per_chr", mode: 'copy', overwrite: true, pattern: "*.sscore"
  publishDir "${params.outdir}/per_chr", mode: 'copy', overwrite: true, pattern: "chr*.scoring.log"

  input:
  tuple val(chrom),
        path(pgen),
        path(pvar),
        path(psam),
        path(weights_list)

  output:
  path("*.chr${chrom}.scored.sscore"), emit: scored_sscores
  path("chr${chrom}.scoring.log"), emit: log

  shell:
  '''
  set -euo pipefail

  exec > >(tee -a chr!{chrom}.scoring.log) 2>&1

  pfile_prefix="pfile_chr!{chrom}"

  echo "=========================================="
  echo "PRS Scoring: ALL TRAITS chr!{chrom}"
  echo "=========================================="
  echo "Start time: $(date)"
  echo "Chrom: !{chrom}"
  echo "Threads: !{task.cpus}"
  echo "weights_list file: !{weights_list}"
  echo "weights.list lines: $(wc -l < !{weights_list})"
  echo "First 10 weights:"
  head -n 10 !{weights_list} || true
  echo ""

  # Link genotype files once
  ln -sf !{pgen} "${pfile_prefix}.pgen"
  ln -sf !{pvar} "${pfile_prefix}.pvar"
  ln -sf !{psam} "${pfile_prefix}.psam"

  if [ ! -s "${pfile_prefix}.pvar" ]; then
    echo "ERROR: PVAR missing/empty after staging: ${pfile_prefix}.pvar"
    exit 1
  fi

  trait_from_weight() {
    local wf="$1"
    local name
    name="$(basename "$wf")"

    for suf in ".tsv.gz" ".csv.gz" ".tsv" ".csv" ".gz"; do
      if [[ "$name" == *"$suf" ]]; then
        name="${name%"$suf"}"
        break
      fi
    done

    name="${name%_clumped}"
    echo "$name"
  }

  # Loop line-by-line through weights.list (robust)
  while IFS= read -r wf; do
    [[ -z "$wf" ]] && continue

    trait="$(trait_from_weight "$wf")"

    echo ""
    echo "------------------------------------------"
    echo "Trait: $trait"
    echo "Weight file: $wf"
    echo "------------------------------------------"

    python3 !{projectDir}/bin/make_scorefile.py \
      --weights "$wf" \
      --chrom !{chrom} \
      --pvar "${pfile_prefix}.pvar" \
      --out_score "${trait}.chr!{chrom}.score.tsv" \
      !{ params.get('allow_strand_flip', false) ? '--allow_strand_flip' : '' }

    if [ ! -s "${trait}.chr!{chrom}.score.tsv" ]; then
      echo "WARNING: Empty score file (no overlaps) for $trait chr!{chrom}"
      touch "${trait}.chr!{chrom}.scored.sscore"
      continue
    fi

    VARIANT_COUNT=$(tail -n +2 "${trait}.chr!{chrom}.score.tsv" | wc -l || true)
    echo "Variants to score: ${VARIANT_COUNT}"

    if [ "${VARIANT_COUNT}" -eq 0 ]; then
      touch "${trait}.chr!{chrom}.scored.sscore"
      continue
    fi

    !{params.get('plink2_bin','plink2')} \
      --pfile "${pfile_prefix}" \
      --score "${trait}.chr!{chrom}.score.tsv" header-read 1 2 3 cols=+scoresums \
      --threads !{task.cpus} \
      --out "${trait}.chr!{chrom}.scored" \
      !{ params.get('write_list_variants', false) ? 'list-variants' : '' }

    if [ ! -f "${trait}.chr!{chrom}.scored.sscore" ]; then
      echo "ERROR: Missing .sscore for $trait chr!{chrom}"
      touch "${trait}.chr!{chrom}.scored.sscore"
    fi

  done < !{weights_list}

  echo ""
  echo "All traits finished for chr!{chrom}."
  echo "End time: $(date)"
  '''
}
