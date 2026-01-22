#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCORE_TRAIT_CHR } from './modules/score_trait_chr'
include { MERGE_TRAIT      } from './modules/merge_trait'
include { BUILD_MATRIX     } from './modules/build_matrix'

/**
 * Trait naming rules:
 * - strip .tsv.gz / .csv.gz / .tsv / .csv / .gz
 * - remove suffix "_clumped"
 */
def traitFromFile(Path f) {
  def name = f.getName()
  for (suf in ['.tsv.gz','.csv.gz','.tsv','.csv','.gz']) {
    if (name.endsWith(suf)) {
      name = name[0..<(name.size()-suf.size())]
      break
    }
  }
  if (name.endsWith('_clumped'))
    name = name[0..<(name.size()-'_clumped'.size())]
  return name
}

workflow {

  // Defensive defaults (prevents null crashes if params not loaded)
  def chrom_list = (params.chroms ?: (1..22)).toList().collect { it.toString() }

  /*
   * 1) WEIGHTS CHANNEL
   * weights_ch emits: tuple(trait, weight_file)
   */
  def weights_ch = Channel
    .fromPath(params.weights_glob, checkIfExists: true)
    .map { f -> tuple(traitFromFile(f), f) }

  /*
   * 2) CHROMS CHANNEL
   * chroms_ch emits: "1".."22"
   */
  def chroms_ch = Channel
    .fromList(chrom_list)

  /*
   * 3) PFILES CHANNEL (staging)
   * pfiles_ch emits: tuple(chrom, pgen, pvar, psam)
   *
   * CRITICAL: PLINK2 cannot read gs:// directly.
   * By declaring them as `file(...)` inputs here, Nextflow stages them locally.
   */
  def pfiles_ch = chroms_ch.map { c ->
    def prefix = "${params.pfile_dir}/${params.pfile_pref}".replace('{chrom}', c)
    tuple(
      c,
      file("${prefix}.pgen"),
      file("${prefix}.pvar"),
      file("${prefix}.psam")
    )
  }

  /*
   * 4) FAN-OUT JOBS: (trait × chrom)
   * trait_chr_jobs emits: tuple(trait, weight_file, chrom)
   */
  def trait_chr_jobs = weights_ch
    .combine(chroms_ch)
    .map { t, w, c -> tuple(t, w, c) }

  /*
   * 5) JOIN JOBS WITH PFILES BY CHROM
   * jobs_with_pfiles emits:
   *   tuple(trait, weight_file, chrom, pgen, pvar, psam)
   */
  def jobs_with_pfiles = trait_chr_jobs
    .map { trait, w, chrom -> tuple(chrom, trait, w) }
    .join(pfiles_ch)   // join on chrom key
    .map { chrom, trait, w, pgen, pvar, psam ->
      tuple(trait, w, chrom, pgen, pvar, psam)
    }

  /*
   * 6) SCORE: per (trait × chrom)
   *
   * IMPORTANT:
   * SCORE_TRAIT_CHR has MULTIPLE emitted outputs (scored_sscore + log).
   * So we must select the emitted channel explicitly.
   */
  def score_res = SCORE_TRAIT_CHR(jobs_with_pfiles)

  // Use ONLY tuple(trait, sscore) output for downstream grouping/merge
  def scored_ch = score_res.scored_sscore

  // (Optional) logs channel if you want to inspect/collect it
  def logs_ch   = score_res.log

  /*
   * 7) GROUP PER TRAIT -> MERGE ACROSS CHROMOSOMES
   * merged_inputs emits: tuple(trait, [sscore1, sscore2, ...])
   */
  def merged_inputs = scored_ch.groupTuple(by: 0)

  /*
   * 8) MERGE: per trait genome-wide table
   * merged emits: path("<trait>.genomewide.sscore.tsv")
   */
  def merged = MERGE_TRAIT(merged_inputs)

  /*
   * 9) BUILD MATRIX ONCE (all traits)
   * collect() makes a single list of all merged trait tables.
   */
  def all_trait_tables = merged.collect()

  BUILD_MATRIX(all_trait_tables)
  
 }

workflow.onComplete {
    log.info """
    Pipeline complete: ${workflow.success}
    workDir : ${workflow.workDir}
    outdir  : ${params.outdir}
    duration: ${workflow.duration}
    """
}

