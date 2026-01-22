nextflow.enable.dsl = 2

process BUILD_MATRIX {
  tag "wide_matrix"

  publishDir "${params.outdir}/prs_merged", mode: 'copy', overwrite: true

  input:
  path(trait_tables)

  output:
  path("prs_matrix_ALL_TRAITS.tsv")

  script:
  """
  set -euo pipefail

  python3 ${projectDir}/bin/build_prs_matrix.py \
    --out prs_matrix_ALL_TRAITS.tsv \
    ${trait_tables.collect{ "--in ${it}" }.join(" ")}

  ls -lh prs_matrix_ALL_TRAITS.tsv
  """
}
