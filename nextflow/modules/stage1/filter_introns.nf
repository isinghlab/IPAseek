process FILTER_INTRONS {
    tag 'filter_introns'
    publishDir "${params.outdir}/stage1/filter_introns", mode: 'copy'
    conda "${projectDir}/../../environment.yml"

    input:
    path annotation_dir

    output:
    path "rnhg38_filtered_introns_cds.rds", emit: intron_annotation

    script:
    """
    Rscript ${projectDir}/../../1_intron_preprocessing/3_filtering_gobj/scripts/1_filtering_gobj.r \\
        --annotation_dir ${annotation_dir} \\
        --outdir ./
    """
}
