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
  # We will build score files for all traits first, then run a single PLINK2 call
  # with --score-list. This avoids reloading the chromosome + recalculating allele
  # frequencies for every single trait.

  : > scorefiles.list
  : > expected_traits.list

  while IFS= read -r wf; do
    [[ -z "$wf" ]] && continue

    trait="$(trait_from_weight "$wf")"

    echo ""
    echo "------------------------------------------"
    echo "Trait: $trait"
    echo "Weight file: $wf"
    echo "------------------------------------------"

    echo "$trait" >> expected_traits.list

    python3 !{projectDir}/bin/make_scorefile.py \
      --weights "$wf" \
      --chrom !{chrom} \
      --pvar "${pfile_prefix}.pvar" \
      --out_score "${trait}.chr!{chrom}.score.tsv" \
      --score_name "$trait" \
      !{ params.get('allow_strand_flip', false) ? '--allow_strand_flip' : '' }

    if [ ! -s "${trait}.chr!{chrom}.score.tsv" ]; then
      echo "WARNING: Empty score file (no overlaps) for $trait chr!{chrom}"
      continue
    fi

    VARIANT_COUNT=$(tail -n +2 "${trait}.chr!{chrom}.score.tsv" | wc -l || true)
    echo "Variants to score: ${VARIANT_COUNT}"

    if [ "${VARIANT_COUNT}" -eq 0 ]; then
      continue
    fi

    # Non-empty score file -> add to --score-list input.
    echo "${trait}.chr!{chrom}.score.tsv" >> scorefiles.list

  done < !{weights_list}

  if [ ! -s scorefiles.list ]; then
    echo "WARNING: No non-empty score files for chr!{chrom}; emitting empty .sscore files"
    while IFS= read -r t; do
      [[ -z "$t" ]] && continue
      touch "${t}.chr!{chrom}.scored.sscore"
    done < expected_traits.list
    exit 0
  fi

  echo ""
  echo "Running ONE PLINK2 scoring call for chr!{chrom} with $(wc -l < scorefiles.list) traits"

  !{params.get('plink2_bin','plink2')} \
    --pfile "${pfile_prefix}" \
    --score-list scorefiles.list header-read 1 2 3 cols=+scoresums \
    --threads !{task.cpus} \
    --out "chr!{chrom}.alltraits.scored" \
    !{ params.get('write_list_variants', false) ? 'list-variants' : '' }

  # Split chrX.alltraits.scored.sscore -> per-trait files with SCORE_SUM column
  python3 !{projectDir}/bin/split_multisscore.py \
    --sscore "chr!{chrom}.alltraits.scored.sscore" \
    --chrom "!{chrom}" \
    --expected_traits expected_traits.list

  echo ""
  echo "All traits finished for chr!{chrom}."
  echo "End time: $(date)"
  '''
}
