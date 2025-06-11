source("/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/rnaseq_pipeline.R")
##############################
# Parse aligner and commands
args <- commandArgs ()
# Set up alignment parameters
N.args <- length(args)
genome.name <- args[(N.args - 4)]
sample.name <- args[(N.args - 3)]
bam.file.location <- args[(N.args - 2)]
save.location <- args[(N.args - 1)]
paired <- args[N.args]

print(genome.name)
print(sample.name)
print(bam.file.location)
print(save.location)
print(paired)

##Get the genome cds object and extract the annotation in relevant format
genome.cds.name <- sprintf("%s.cds.rds", genome.name)
annotation.dir <- ("/scratch/user/richa.rashmi.1202/genomes/annotation_objects")
genome.cds.location <- file.path(annotation.dir, genome.cds.name)
ls.annotations <- ExtractAnnotation(genome.cds.location, genome.name, num.cores = 12)
grl.gene.annotation <- ls.annotations$gene.annotation
gr.exon.annotation <- ls.annotations$exon.annotation

ls.counts <- ReadRNASeqBam(genome = genome.name, bam.file.location, grl.gene.annotation, 
	gr.exon.annotation, paired = paired)

m.gene.counts <- as.matrix(ls.counts$gene.counts)
m.exon.counts <- as.matrix(ls.counts$exon.counts)
lib.size <- ls.counts$library.size

df.lib.size <- DataFrame(library.size = lib.size)
rownames(df.lib.size) <- sample.name
colnames(m.gene.counts) <- sample.name
colnames(m.exon.counts) <- sample.name

se.gene <- SummarizedExperiment(rowData = grl.gene.annotation, colData = df.lib.size,
  assays = SimpleList(count = m.gene.counts))
se.exon <- SummarizedExperiment(rowData = gr.exon.annotation, colData = df.lib.size,
  assays = SimpleList(count = m.exon.counts))

saveRDS(se.gene, file = file.path(save.location, sprintf("%s_gene_counts.rds", sample.name)))
saveRDS(se.exon, file = file.path(save.location, sprintf("%s_exon_counts.rds", sample.name)))
