#!/bin/bash
#SBATCH --job-name=starindex         # Job name
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=8          #Request 8 tasks/cores per node
#SBATCH --mem=56G                     #Request 56GB per node 
#SBATCH --time=24:00:00              #Time limit hrs:min:sec
#SBATCH --output=starindexgenerate_%j.log   # Standard output and error log
pwd; hostname; date

# Launch my program.
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch/user/isingh/genomes/Homo_sapiens/UCSC/hg38/Sequence/STARIndex --genomeFastaFiles /scratch/user/isingh/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /scratch/user/isingh/genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf --sjdbGTFfeatureExon exon --sjdbOverhang 100
