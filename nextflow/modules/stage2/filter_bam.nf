process FILTER_BAM {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage2/filter_bam/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(bam), path(bai)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_uniq.bam"), path("${sample_id}_uniq.bam.bai"), emit: uniq_bam

    script:
    """
    samtools view -h -@ ${task.cpus} ${bam} \\
        | grep -E '^@|NH:i:1\\b' \\
        | samtools view -bS -@ ${task.cpus} - \\
        > ${sample_id}_uniq.bam

    samtools index -@ ${task.cpus} ${sample_id}_uniq.bam
    """
}
