#!/usr/bin/env nextflow


process t1 {

	publishDir path: '/home/yiyuan/Data/Nextflow/test', mode:"copy", overwrite: true, saveAs: {

		switch (it) {
			case { it.endsWith(".txt")}:
				return 'a.txt'
			case { it.endsWith(".csv")}:
				return it
			default:
				println "not matched: " + it; exit 1;
		}

	}


	echo true

	output:
	file 'file.txt' into a
	file 'file.csv' into b

	"""
		echo 'txt' > file.txt
		echo 'csv' > file.csv
	"""
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
