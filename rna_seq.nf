#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process trimFastp {
    
    label 'fastp'
    publishDir "${params.output}/trim", mode: 'copy'
    container "biocontainers/fastp:v0.20.1_cv1"

    input:
        tuple val(runID), path(fq1), path(fq2), val(sampleID)
    output:
        tuple val(sampleID), path("${sampleID}.trim.R1.fq.gz"), path("${sampleID}.trim.R2.fq.gz"), emit: trimReads
        path("${sampleID}.trim_report.*"), emit: trimReports
    script:
    """
    fastp -i ${fq1} -I ${fq2} -o ${sampleID}.trim.R1.fq.gz -O ${sampleID}.trim.R2.fq.gz -w ${params.threads} -h ${sampleID}.trim_report.html -j ${sampleID}.trim_report.json -q 25 -l 30 --length_limit 0 -g -x
    """

}

process alignHisat {

    label 'hisat2'
    publishDir "${params.output}/align", mode: 'copy'
    container "nanozoo/hisat2:2.2.1.commit7e01700--5e923e8"

    maxForks 1

    input:
        tuple val(sampleID), path(trim1), path(trim2), path(hisatIdx)
    output:
        tuple val(sampleID), path("${sampleID}.bam"), emit: alignBam
        path("${sampleID}.alignerout.txt"), emit: alignReport
    script:
    """
    tar -xzvf ${hisatIdx}

    hisat2 -p ${params.threads} --rg-id ${sampleID} --rg LB:tx --rg PL:illumina --rg PU:barcode --rg SM:${sampleID} --dta -x hisat2_index/genome -1 ${trim1} -2 ${trim2} -S out.sam --summary-file ${sampleID}.alignerout.txt

    samtools view -1 --threads ${params.threads} -o out.bam out.sam

    samtools sort -@ ${params.threads} -O BAM -o ${sampleID}.bam out.bam

    samtools index -@ ${params.threads} ${sampleID}.bam
    """

}

process featureCounts {

    label 'featureCounts'
    publishDir "${params.output}/counts", mode: 'copy'
    container "pegi3s/feature-counts:2.0.0"

    input:
        tuple val(sampleID), path(bam), path(annotation)
    output:
        path("${sampleID}.cts.txt"), emit: countTxt
        path("${sampleID}.cts.txt.*"), emit: countStats
    script:
    """
    gunzip -f ${annotation}

    featureCounts -s 0 -M --fraction -J --ignoreDup -T ${params.threads} -p -g gene_id -a gencode.v48.annotation.gtf -o ${sampleID}.cts.txt ${bam}
    """

}

process transcriptCounts {

    label 'stringtie'
    publishDir "${params.output}/counts", mode: 'copy'
    container "mgibio/stringtie:v2.2.1-focal"

    input:
        tuple val(sampleID), path(bam), path(transcript)
    output:
        path("${sampleID}_stringtie.tar.gz"), emit: stringtieCounts
    script:
    """
    gunzip -f ${transcript}

    mkdir -p ${sampleID}_stringtie

    stringtie ${bam} -p ${params.threads} -G transcripts.gtf -B -e -o ${sampleID}.denovo.gtf -A ${sampleID}.fpkm.txt

    mv *.ctab ${sampleID}_stringtie

    mv ${sampleID}.denovo.gtf ${sampleID}_stringtie

    mv ${sampleID}.fpkm.txt ${sampleID}_stringtie

    tar -czvf ${sampleID}_stringtie.tar.gz ${sampleID}_stringtie
    """

}

process statAnalysis {

    label 'deseq'
    publishDir "${params.output}/figures", mode: 'copy'

    input:
        path(countTables)
        path(metadata)
        path(deq)
    output:
        tuple path("volcano.png"), path("pca.png"), path("heatmap.png"), path("distance_heatmap.png"), path("GO_analysis.csv") emit: figures
    script:
    """
    Rscript ${deq}
    """

}

workflow {

    // Input raw fastq data
    Channel.fromFilePairs("${params.inputFQ}/*_{1,2}.fastq.gz", flat: true)
    .set{ runIDMap }

    // Input design file
    Channel.fromPath("${params.designFile}")
    .splitCsv( header: true, sep: "\t" )
    .map { row -> tuple( row.RunID, row.SampleID ) }
    .set{ sampleIDMap }

    // Map Sample ID to Run ID
    runIDMap.join( sampleIDMap )
    .set{ trimInput }
    
    // Trim
    trimFastp(trimInput)

    // Align
    alignHisat(trimFastp.out.trimReads.combine(Channel.fromPath(params.hisatIdx)))

    // Feature Counts
    featureCounts(alignHisat.out.alignBam.combine(Channel.fromPath(params.annotationFile)))

    // Transcript Counts
    transcriptCounts(alignHisat.out.alignBam.combine(Channel.fromPath(params.transcriptFile)))

    // Stat Analysis
    Channel.fromPath("${params.designFile}")
    .set{ metadata }

    Channel.fromPath("${params.statAnalysisFile}")
    .set{ deq }

    statAnalysis(featureCounts.out.countTxt.collect(), metadata, deq)

}