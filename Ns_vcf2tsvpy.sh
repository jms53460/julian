#!/bin/bash
#SBATCH --job-name=Ns_vcf2tsvpy                         # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=6:00:00                                  # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_vcf2tsvpy.out        # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_vcf2tsvpy.err         # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460
ml Miniconda3/23.5.2-0
source activate /home/jms53460/vcf2tsvpy
vcf2tsvpy --input_vcf Ns.vcf.gz --out_tsv Ns.vcf_table.tsv 
conda deactivate