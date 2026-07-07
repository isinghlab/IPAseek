process STAR_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/stage2/star_align/${sample_id}", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path star_index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${sample_id}_

    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}.bam
    """
}
