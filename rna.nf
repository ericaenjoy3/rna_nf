#!/usr/bin/env nextflow

import java.nio.file.*

/*
 * 'rna' was implemented by Yiyuan Liu of mssm.
 */

/*
 * 'rna': A nextflow-based RNA-seq analysis pipeline from RNA sequencing data
 */

/* requirement:
 * - STAR/tophat2/bowtie2
 * - samtools/sambamba
 */

/*
 * Nextflow Version check
 */
try {
    if( !nextflow.version.matches('0.26+') ){
        throw GroovyException('This workflow requires Nextflow version 0.26 or greater.')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version 0.26 or greater required! You are running v$nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * Miscellaneous code for the pipeline
 */
boolean isSingleFile(object) {
    object instanceof Path
}

String getDir(String dir, String subdir = 'base_dir', params, canonPath = false) {
	if (canonPath) {
		return new File(params[dir].get(subdir, subdir)).getCanonicalPath()
	} else {
		return params[dir].get(subdir, subdir)
	}
}

String concatDir(String subdir1, String subdir2) {
  new File(subdir1, subdir2).getPath()
}

String joinPath(String first, String... args){
	Path path = Paths.get(first, args)
	path.toString()
}

def action = { params, l1, l2, indicator ->
	switch (indicator){
		case ~/(string|integer)/:
			return params[l1].get(l2, '')
		case ~/file/:
			return new File(params[l1].get(l2, '')).getCanonicalPath()
		case ~/boolean/:
			return params[l1].get(l2, 'false').toBoolean()
		case ~/(list|map)/:
			return params[l1].get(l2, '')
		default:
			throw new IllegalArgumentException( "Can't accept the indicator value: " + indicator)
	}
}

def testBoolean = { it -> "true".equals(it) || "false".equals(it) }
def testInteger = { it -> it.toString().isInteger()}

/*
 * Validate setting parameters
 */
// leaf keys carry data types and values are verified against specified data type.
def validSetting(String k, Object v) {

	(full, type) = (k =~ /.*_+([^_]+)\\/?$/)[0]

	println "parameter: " + k + ", required type: " + type

	switch (type) {

		case 'file':

			if (!(new File(v).exists())) {
				throw new IllegalArgumentException( "File does not exist: " + v)
			}
			break

		case {it == 'dir' && k =~ /base_dir/}:

			if (!(new File(v).getParentFile().exists())) {
				throw new IllegalArgumentException("Directory does not exists: " + new File(v).getParent())
			}
			break

		case {it == 'dir' && !(k =~ /base_dir/)}:

			if (!(v)) {
				throw new IllegalArgumentException( "Parameters unset: " + k)
			}
			break

		case 'string':

			if (!(v instanceof String)) {
				throw new IllegalArgumentException( "Parameters should be string: " + k)
			}
			break

		case 'integer':

			if (!(v instanceof Integer)) {
				throw new IllegalArgumentException( "Parameters should be integer: " + k)
			}
			break

		case 'boolean':

			if (!(v instanceof Boolean)) {
				throw new IllegalArgumentException( "Parameters should be boolean: " + k)
			}
			break

		default:
			throw new IllegalArgumentException( "Can't handle key: " + k + "; and value: " + v)
	}
}

def call_validSetting(m) {

   m.findResults { k, v -> v instanceof Map ? call_validSetting(v) : validSetting(k, v) }

}

params.keySet().findAll { it -> !(it =~ /(?i)input/) }.each {k -> call_validSetting(params[k])}

/*
 * sanity check for samtools_sambamba_string
 */
if (!(action.call(params, 'tools', 'samtools_sambamba_string', 'string') in ['sambamba', 'samtools'])) {

	throw new IllegalArgumentException( "samtools_sambamba_string parameter can only be samtools or sambamba: " + action.call(params, 'tools', 'samtools_sambamba_string', 'string') + " not accepted.")

}

/*
 * Process reads
 */
LIB_RUN_LOCAL_FASTQS = Channel.from(
  action.call(params, 'input', 'raw_fastq_dir', 'map').collect {k, v -> [k,v]})

LIB_RUN_LOCAL_FASTQS
  .map { v -> [v[0], v[1].collect{file(it)}] }
  .into {SALMON_FQ; STAR_FQ; TOPHAT_FQ, TEST_FQ}

TEST_FQ.subscribe {println it}

TRANS_FA = Channel.fromPath(action.call(params, 'fasta', 'transcript', 'file'))
	.ifEmpty { exit 1, "Index not found: " + action.call(params, 'fasta', 'transcript', 'file')}
(SALMON_TRANS_FA) = TRANS_FA.into(1)

GENOME_FA = Channel.fromPath(action.call(params, 'fasta', 'genome', 'file'))
	.ifEmpty { exit 1, "Index not found: " + action.call(params, 'fasta', 'genome', 'file')}
(STAR_GENOME_FA, BOWTIE2_GENOME_FA,
	TOPHAT_GENOME_FA) = GENOME_FA.into(3)

GTF = Channel.fromPath(action.call(params, 'annotation', 'gtf', 'file'))
	.ifEmpty { exit 1, "Index not found: " + action.call(params, 'annotation', 'gtf', 'file')}
(STAR_GTF, TOPHAT_GTF) = GTF.into(2)

if (action.call(params, 'tools', 'samtools_sambamba_string', 'string')) {
	bamIndex_tool = Channel.from(action.call(params, 'tools', 'samtools_sambamba_string', 'string'))
	(STAR_BAM_INDEX_TOOL, TOPHAT_BAM_INDEX_TOOL) = bamIndex_tool.into(2)
}
/*
 * Step 1: Build STAR/tophat/salmon/kallisto index if non-existent
 */
// salmon index process
process salmonIndex {

	label 'multiCore'

  tag "salmonIndex"

  publishDir path: getDir('reference', 'base_dir', params, canonPath = true), mode:"copy", saveAs: {
		return getDir('reference', 'salmon_dir', params, canonPath = false)
	}

  input:
  file fa from SALMON_TRANS_FA

  output:
  file "salmon_idx" into SALMON_INDEX

	when:
	action.call(params, 'compute', 'run_salmon_index_boolean', 'boolean')

  script:
	gencode = action.call(params, 'annotation', 'gencode_boolean', 'boolean')

	if (!testBoolean.call(gencode)) {
		throw new IllegalArgumentException("gencode_boolean (" + gencode + ") is not Boolean")
	}

  switch (gencode) {
    case true:
      gencode_flag = '--gencode'
      break
    default:
      gencode_flag = ''
      break
  }

  threads = task.cpus

  """
  salmon index -t ${fa} -i salmon_idx --type quasi -k 31 -p ${threads} ${gencode_flag}
  """
}

// star index process
process startIndex {

	label 'multiCore'

	tag "starIndex"

  publishDir path: getDir('reference', 'base_dir', params, canonPath = true), mode:"copy", saveAs: {
		return getDir('reference', 'star_dir', params, canonPath = false)
	}

  input:
  file fa from STAR_GENOME_FA
	file gtf from STAR_GTF

  output:
  file "star_idx" into STAR_INDEX

	when:
	action.call(params, 'compute', 'run_star_index_boolean', 'boolean')

  script:
	sjdbOverhang = action.call(params, 'reference', 'star_sjdbOverhang_integer', 'integer')

	if (!testInteger.call(sjdbOverhang)) {
		throw new IllegalArgumentException("sjdbOverhang ("sjdbOverhang + ") is not Integer")
	}

  threads = 1
	"""
		mkdir star_idx
		STAR \
				--runMode genomeGenerate \
				--runThreadN ${threads} \
				--sjdbGTFfile ${gtf} \
				--sjdbOverhang ${sjdbOverhang} \
				--genomeDir star_idx \
				--genomeFastaFiles ${fa}
	"""
}

// bowtie2 index process
process bowtie2Index {

	label 'multiCore'

	tag "bowtie2Index"

	publishDir path: getDir('reference', 'base_dir', params, canonPath = true), mode:"copy", saveAs: {
		return getDir('reference', 'tophat_bowtie2_dir', params, canonPath = false)
	}

  input:
  file fa from BOWTIE2_GENOME_FA

  output:
  file "genome_bt2.*" into BOWTIE2_INDEX

	when:
	action.call(params, 'compute', 'run_bowtie2_index_boolean', 'boolean')

  script:
	threads = task.cpus
  """
  	bowtie2-build --threads ${threads} ${fa} genome_bt2
  """
}

/*
* Step 2: alignment of raw reads
*/
// starAlign process
if (!action.call(params, 'compute', 'run_star_index_boolean', 'boolean') &&
	action.call(params, 'compute', 'run_star_boolean', 'boolean')) {

		println "star index was built already"

		def idx_filename = concatDir(getDir('reference', 'base_dir', params, canonPath = true),
		getDir('reference', 'star_dir', params, canonPath = false))

		println "pre-built star index: " + idx_filename

		STAR_INDEX = Channel.fromPath("${idx_filename}")
			.ifEmpty { exit 1, "Index not found: ${idx_filename}" }
}

process starAlign {

	label 'multiCore'

	tag "starAlign: $sample"

	publishDir path: joinPath(getDir('output', 'base_dir', params, canonPath = true),
		getDir('output', 'star_dir', params, canonPath = false)),
		mode: "copy", saveAs: {

			switch (it) {
				case { it.endsWith('.bam' ) }:
					return joinPath("$sample", "${sample}.bam"
				case { it.endsWith('.bam.bai')}:
					return joinPath("$sample", "${sample}.bam.bai")
				default:
					return "$sample" + it
			}

		}

  input:
  file index from STAR_INDEX.collect()
  set val(sample), file(reads) from STAR_FQ
	val tool from STAR_BAM_INDEX_TOOL

  output:
  set val(sample), file("${sample}_Aligned.sortedByCoord.out.bam") into STAR_BAM

	when:
	action.call(params, 'compute', 'run_star_boolean', 'boolean')

  script:
  def single = reads instanceof Path
  threads = task.cpus

	switch (tool){
		case 'sambamba':
			flag = '--nthreads ' + threads
			break
		case 'samtools':
			flag = ''
			break
		default:
			throw new IllegalArgumentException( "Can't accept the samtools_sambamba_string: " + toools)
	}

  if (single) {
    """

			echo "\n===\nProcessing starAlign\n===\nread: ${reads}\n===\n"

			STAR --runThreadN ${threads} \
				 --twopassMode Basic \
				 --genomeDir ${index} \
				 --readFilesIn ${reads} \
				 --readFilesCommand zcat \
				 --outSAMtype BAM SortedByCoordinate \
				 --chimSegmentMin 20 \
				 --outFilterIntronMotifs RemoveNoncanonical \
				 --outFilterMultimapNmax 20 \
				 --alignIntronMin 20 \
				 --alignIntronMax 1000000 \
				 --alignMatesGapMax 1000000 \
				 --outFilterType BySJout \
				 --alignSJoverhangMin 8 \
				 --alignSJDBoverhangMin 1 \
				 --outFileNamePrefix "${sample}_"

			${tool} index ${flag} ${sample}_Aligned.sortedByCoord.out.bam
    """
  } else {
    """

	    echo "\n===\nProcessing starAlign\n===\nR1: ${reads[0]} and R2: ${reads[1]}\n===\n"

			STAR --runThreadN ${threads}  \
					 --twopassMode Basic --genomeDir ${index} \
					 --readFilesIn ${reads[0]} ${reads[1]} \
					 --readFilesCommand zcat \
					 --outSAMtype BAM SortedByCoordinate \
					 --chimSegmentMin 20 \
					 --outFilterIntronMotifs RemoveNoncanonical \
					 --outFilterMultimapNmax 20 \
					 --alignIntronMin 20 \
					 --alignIntronMax 1000000 \
					 --alignMatesGapMax 1000000 \
					 --outFilterType BySJout \
					 --alignSJoverhangMin 8 \
					 --alignSJDBoverhangMin 1 \
					 --outFileNamePrefix "${sample}_"

			${tool} index ${flag} ${sample}_Aligned.sortedByCoord.out.bam
    """
  }

}

// tophatAlign process
if (!action.call(params, 'compute', 'run_tophat_bowtie2_index_boolean', 'boolean') &&
	action.call(params, 'compute', 'run_tophat_boolean', 'boolean')) {

		println "bowtie2 index was built already"

		def idx_filename = concatDir(getDir('reference', 'base_dir', params, canonPath = true),
		getDir('reference', 'tophat_bowtie2_dir', params, canonPath = false))

		println "pre-built bowtie2 index: " + idx_filename

		BOWTIE2_INDEX = Channel.fromPath("${idx_filename}")
			.ifEmpty { exit 1, "Index not found: ${idx_filename}" }
}

process tophatAlign {

	label 'multiCore'

	tag "tophatAlign: $sample"

	publishDir path: joinPath(getDir('output', 'base_dir', params, canonPath = true),
		getDir('output', 'tophat_dir', params, canonPath = false)),
		mode: "copy", saveAs: {

			switch (it) {
				case { it.endsWith('.bam' ) }:
					return joinPath("$sample", "${sample}.bam"
				case { it.endsWith('.bam.bai')}:
					return joinPath("$sample", "${sample}.bam.bai")
				default:
					return "$sample" + it
			}

		}

  input:
  file index from BOWTIE2_INDEX.collect()
	file fa from TOPHAT_GENOME_FA
	file gtf from TOPHAT_GTF
  set val(sample), file(reads) from TOPHAT_FQ
	val tool from TOPHAT_BAM_INDEX_TOOL

  output:
  set val(sample), file("accepted_hits.bam") into TOPHAT_BAM

	when:
	action.call(params, 'compute', 'run_tophat_boolean', 'boolean')

  script:
	threads = task.cpus
	strand = action.call(params, 'input', 'tophat_strand', 'string')

	switch (tool){
		case 'sambamba':
			flag = '--nthreads ' + threads
			break
		case 'samtools':
			flag = ''
			break
		default:
			throw new IllegalArgumentException( "Can't accept the samtools_sambamba_string: " + toools)
	}

	if (!(strand instanceof String && strand in ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand'])) {
		throw new IllegalArgumentException("tophat_strand (" +
		strand + ") is not in fr-unstranded, fr-firststrand or fr-secondstrand")
	}

	index_base = index[0].toString() - ~/.\d.bt2/
  def single = reads instanceof Path
  threads = task.cpus
  if (single) {
    """

			echo "\n===\nProcessing tophatAlign\n===\nread: ${reads}\n===\n"

			tophat -p ${threads} \
				-G ${gtf} \
				--no-novel-juncs \
				-o ${sample} \
				--library-type ${strand} \
				${index_base} ${reads}

			${tool} index ${flag} accepted_hits.bam
    """
  } else {
    """

	    echo "\n===\nProcessing tophatAlign\n===\nR1: ${reads[0]} and R2: ${reads[1]}\n===\n"

			tophat -p ${threads} \
				-G ${gtf} \
				--no-novel-juncs \
				-o ${sample} \
				--library-type ${strand} \
				${index_base} ${reads[0]} ${reads[1]}

			${tool} index ${flag} accepted_hits.bam
    """
  }
}

/*
* Step 3: bam QC
*/
// collect bam files
STAR_BAM
	.mix(TOPHAT_BAM)
	.set{BAM_NO_QC}

/*
* Step 4: salmon quantification
*/
// salmon quantification process
// run_salmon_index (true): run salmon index
// run salmon_index (false) && run salmon (true): load salmon index
if (!action.call(params, 'compute', 'run_salmon_index_boolean', 'boolean') &&
	action.call(params, 'compute', 'run_salmon_boolean', 'boolean')) {

		println "salmon index was built already"

		def idx_filename = concatDir(getDir('reference', 'base_dir', params, canonPath = true),
		getDir('reference', 'salmon_dir', params, canonPath = false))

		println "pre-built salmon index: " + idx_filename

		SALMON_INDEX = Channel.fromPath("${idx_filename}")
			.ifEmpty { exit 1, "Index not found: ${idx_filename}" }

}

process salmonQuant {

	label 'multiCore'

  tag "salmonQuant: $sample"

	publishDir path: concatDir(getDir('output', 'base_dir', params, canonPath = true),
		getDir('output', 'salmon', params, canonPath = false)),
		mode: "copy"

  input:
  file index from SALMON_INDEX.first()
  set val(sample), file(reads) from SALMON_FQ

  output:
  set val(sample), file("$sample") into SALMON_OUT

	when:
	action.call(params, 'compute', 'run_salmon_boolean', 'boolean')

  script:
  def single = reads instanceof Path
  threads = 1
  if (single) {
    """

	    echo "\n===\nProcessing salmonQuant\n===\nread: ${reads}\n===\n"

	    salmon quant --quiet -i ${index} -l A -r ${reads} -o "$salmon" -p ${threads}
    """
  } else {
    """

	    echo "\n===\nProcessing salmonQuant\n===\nR1: ${reads[0]} and R2: ${reads[1]}\n===\n"

	    salmon quant --quiet -i ${index} -l A -1 ${reads[0]} -2 ${reads[1]} -o "$sample" -p ${threads}
    """
  }
}

workflow.onComplete {

	log.info "Pipeline execution summary"
	log.info "---------------------------"
	log.info "Completed at: ${workflow.complete}"
	log.info "Duration    : ${workflow.duration}"
	log.info "Success     : ${workflow.success}"
	log.info "workDir     : ${workflow.workDir}"
	log.info "exit status : ${workflow.success ? 'OK' : 'failed' }"
	log.info "Error report: ${workflow.errorReport ?: '-'}"

}

workflow.onError {
	log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
