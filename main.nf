#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCORE_CHR_ALL_TRAITS } from './modules/score_chr_all_traits'
include { MERGE_TRAIT         } from './modules/merge_trait'
include { BUILD_MATRIX        } from './modules/build_matrix'


def parseChroms(def chromsParam) {
  if (chromsParam == null)
    return (1..22).collect { it.toString() }

  if (chromsParam instanceof List)
    return chromsParam.collect { it.toString() }

  def s = chromsParam.toString().trim()

  if (s.contains('..')) {
    def parts = s.split(/\.\./)
    int a = parts[0].toInteger()
    int b = parts[1].toInteger()
    return (a..b).collect { it.toString() }
  }

  if (s.contains(',')) {
    return s.split(',').collect { it.trim() }.findAll { it }
  }

  return [ s ]
}

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

def traitFromSscore(Path f) {
  def name = f.getName()
  name = name.replaceFirst(/\.chr\d+\.scored\.sscore$/, '')
  return name
}

workflow {

  def chrom_list   = parseChroms(params.get('chroms', null))
  def pfile_dir    = params.get('pfile_dir',  '../urvashi_analysis/pgen_dir')
  def pfile_pref   = params.get('pfile_pref', 'snpqc.chr{chrom}')
  def weights_glob = params.get('weights_glob', '../urvashi_analysis/inputs/test/*.tsv.gz')

  /*
   * weights_ch emits: tuple(trait, weight_file)
   */
  def weights_ch = Channel
    .fromPath(weights_glob, checkIfExists: true)
    .map { f -> tuple(traitFromFile(f), file(f)) }

  /*
   * Optional trait filter: --traits A,B,C
   */
  if (params.get('traits', null)) {
    def keep = (params.traits instanceof List
      ? params.traits
      : params.traits.toString().split(',').collect { it.trim() }
    ) as Set

    weights_ch = weights_ch.filter { trait, wf -> keep.contains(trait) }
  }

  /*
   * Create ONE weights.list file: one absolute path per line
   */

  def weights_list_ch = weights_ch
    .map { trait, wf -> wf.toAbsolutePath().toString() }
    .unique()
    .collectFile(name: 'weights.list', newLine: true)


  /*
   * Chroms channel
   */
  def chroms_ch = Channel.fromList(chrom_list)

  /*
   * Pfiles channel: tuple(chrom, pgen, pvar, psam)
   */
  def pfiles_ch = chroms_ch.map { c ->
    def prefix = "${pfile_dir}/${pfile_pref}".replace('{chrom}', c)
    tuple(
      c,
      file("${prefix}.pgen"),
      file("${prefix}.pvar"),
      file("${prefix}.psam")
    )
  }

  /*
   * One job per chromosome:
   * tuple(chrom, pgen, pvar, psam, weights.list)
   */
  def chr_jobs = pfiles_ch
    .combine(weights_list_ch)
    .map { chrom, pgen, pvar, psam, weights_list ->
      tuple(chrom, pgen, pvar, psam, weights_list)
    }

  /*
   * Score all traits within each chromosome job
   */
  def score_res = SCORE_CHR_ALL_TRAITS(chr_jobs)

  /*
   * Flatten -> group -> merge per trait
   */
  def scored_flat  = score_res.scored_sscores.flatten()
  def scored_pairs = scored_flat.map { f -> tuple(traitFromSscore(f), f) }
  def merged_inputs = scored_pairs.groupTuple(by: 0)

  def merged = MERGE_TRAIT(merged_inputs)

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
