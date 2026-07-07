include { FILTER_INTRONS } from '../modules/stage1/filter_introns'

workflow STAGE1 {
    take:
    annotation_dir  // path to annotation directory

    main:
    FILTER_INTRONS(annotation_dir)

    emit:
    intron_annotation = FILTER_INTRONS.out.intron_annotation
}
