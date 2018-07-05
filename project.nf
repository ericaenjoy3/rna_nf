#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


// the Channel the location of Raw Data fastqs):
LIB_RUN_LOCAL_FASTQS = Channel.from(
    params.input.raw_reads_paths.collect {
        k, v -> v.collect{k2, v2 -> [k,k2]+v2}}.sum())

LIB_RUN_LOCAL_FASTQS.map{ v -> [v[0], v[1], file(v[2]), file(v[3])]}.set{ LIB_RUN_FASTQS }

/*
 * FastQC the input files
 */

LIB_RUN_FASTQS_FOR_QC = Channel.create()
LIB_RUN_FASTQS
    .tap(LIB_RUN_FASTQS_FOR_QC)
    .set{LIB_RUN_FASTQS}

// LIB_RUN_FASTQS_FOR_QC.subscribe { println it }

LIB_RUN_FASTQS_FOR_QC
    .filter { it -> params.get('do_fastqc', 'false').toBoolean() }
    .map{ v -> [v[0], v[1], [[1,file(v[2])], [2,file(v[3])]]]}
    .flatMap{
        vs -> vs[2].collect{
            it -> [vs[0],
                   vs[1],
                   it[0],
                   it[1]] } }
    .set {LIB_RUN_SIDE_FASTQS_FOR_QC}

println "++++++++++++"
LIB_RUN_SIDE_FASTQS_FOR_QC.println()
