nextflow.enable.dsl = 2

process SCORE_TRAIT_CHR {
  tag { "${trait}:chr${chrom}" }

  publishDir "${params.outdir}/per_chr", mode: 'copy', overwrite: true, pattern: "*.sscore"

  input:
  tuple val(trait), path(weight_file), val(chrom), path(pgen), path(pvar), path(psam)

  output:
  tuple val(trait), path("${trait}.chr${chrom}.scored.sscore"), emit: scored_sscore
  path("${trait}.chr${chrom}.scoring.log"), emit: log

  script:
  def pfile_prefix = "pfile_chr${chrom}"

  """
  set -euo pipefail

  exec > >(tee -a ${trait}.chr${chrom}.scoring.log) 2>&1

  echo "=========================================="
  echo "PRS Scoring: ${trait} chr${chrom}"
  echo "=========================================="
  echo "Start time: \$(date)"
  echo "Weight file: ${weight_file}"
  echo "PGEN staged: ${pgen}"
  echo "PVAR staged: ${pvar}"
  echo "PSAM staged: ${psam}"
  echo "PFILE prefix: ${pfile_prefix}"
  echo "Threads: ${task.cpus}"
  echo ""

  ln -sf ${pgen} ${pfile_prefix}.pgen
  ln -sf ${pvar} ${pfile_prefix}.pvar
  ln -sf ${psam} ${pfile_prefix}.psam

  if [ ! -s "${pfile_prefix}.pvar" ]; then
    echo "ERROR: PVAR missing/empty after staging: ${pfile_prefix}.pvar"
    echo "FATAL: empty score file for ${trait} chr${chrom}" >&2
    exit 1
  fi

  echo "Step 1: Creating score file from weights..."
  python3 ${projectDir}/bin/make_scorefile.py \
    --weights ${weight_file} \
    --chrom ${chrom} \
    --pvar ${pfile_prefix}.pvar \
    --out_score ${trait}.chr${chrom}.score.tsv \
    ${params.allow_strand_flip ? '--allow_strand_flip' : ''}

  if [ ! -s ${trait}.chr${chrom}.score.tsv ]; then
    echo "WARNING: Empty score file (no overlaps)."
    touch ${trait}.chr${chrom}.scored.sscore
    exit 0
  fi

  VARIANT_COUNT=\$(tail -n +2 ${trait}.chr${chrom}.score.tsv | wc -l)
  echo "Variants to score: \${VARIANT_COUNT}"

  if [ "\${VARIANT_COUNT}" -eq 0 ]; then
    touch ${trait}.chr${chrom}.scored.sscore
    exit 0
  fi

  echo "Step 2: Running PLINK2 scoring..."
  ${params.plink2_bin} \
    --pfile ${pfile_prefix} \
    --score ${trait}.chr${chrom}.score.tsv header-read 1 2 3 cols=+scoresums \
    --threads ${task.cpus} \
    --out ${trait}.chr${chrom}.scored \
    ${params.write_list_variants ? 'list-variants' : ''}

  if [ -f "${trait}.chr${chrom}.scored.sscore" ]; then
    SAMPLE_COUNT=\$(tail -n +2 ${trait}.chr${chrom}.scored.sscore | wc -l)
    echo "SUCCESS: Scored \${SAMPLE_COUNT} samples"
    ls -lh ${trait}.chr${chrom}.scored.sscore
  else
    echo "ERROR: PLINK2 did not produce sscore file"
    touch ${trait}.chr${chrom}.scored.sscore
  fi

  echo "End time: \$(date)"
  """
}
