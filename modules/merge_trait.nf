nextflow.enable.dsl = 2

process MERGE_TRAIT {
  tag "${trait}"

  publishDir "${params.outdir}/prs_merged", mode: 'copy', overwrite: true

  input:
  tuple val(trait), path(sscore_files)

  output:
  path("${trait}.genomewide.sscore.tsv")

  script:
  """
  set -euo pipefail

  python3 ${projectDir}/bin/merge_sscore_parts.py \
    --trait ${trait} \
    --out ${trait}.genomewide.sscore.tsv \
    ${sscore_files.collect{ "--part ${it}" }.join(" ")}

  ls -lh ${trait}.genomewide.sscore.tsv
  """
}
