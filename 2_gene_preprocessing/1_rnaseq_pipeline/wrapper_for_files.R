# Parse aligner and commands
##Sample Run: Rscript /cbio/cllab/home/singh/clusterRuns/gene_expression_counting/wrapper_for_files.R /cbio/cllab/home/singh/landau2015_data
##Sample Run: Rscript /cbio/cllab/home/singh/clusterRuns/gene_expression_counting/wrapper_for_files.R /cbio/cllab/nobackup/singh/CLL_RNA_Seq
##Sample Run: Rscript /cbio/cllab/home/singh/clusterRuns/gene_expression_counting/wrapper_for_files.R /cbio/cllab/nobackup/singh/mm_rna_seq
##Sample Run: Rscript /cbio/cllab/home/singh/clusterRuns/gene_expression_counting/wrapper_for_files.R /cbio/cllab/nobackup/singh/tamu_data_december2017/17135GAH_N17196

args <- commandArgs ()
# Get the location of the file
bam_at_this_location <- args[length (args)]
print(bam_at_this_location)

all.bam.files <- list.files(path =  bam_at_this_location, pattern = "*.bam$", full.names = TRUE, recursive = T)
for(bam_file in all.bam.files){
	qsub.cmd <- sprintf("qsub -v bam_file=%s /cbio/cllab/home/singh/clusterRuns/gene_expression_counting/run_gene_counts.sh", bam_file)
	system(qsub.cmd)
}
