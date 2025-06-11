#!/bin/bash
#SBATCH --job-name=count         # Job name
#SBATCH --ntasks=1                   # Run a single task	
#SBATCH --cpus-per-task=12           # Number of CPU cores per task
#SBATCH --mem=56G                   # Job memory request
#SBATCH --time=02:00:00              # Time limit hrs:min:sec
#SBATCH --output=2_gene_preprocessing/3_gene_expression/logs/count_%j.log   # Standard output and error log
pwd; hostname; date

conda activate wgbs

# Launch my program.
/scratch/user/richa.rashmi.1202/.conda/envs/wgbs/bin/R --vanilla --file=/scratch/user/richa.rashmi.1202/ipa/ipa_pipeline/2_gene_preprocessing/1_rnaseq_pipeline/countGeneExpression.R ${genome_name} ${sample_name} ${bam_file} ${save_location} ${paired}
