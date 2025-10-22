#!/bin/bash
#SBATCH --job-name=Zm_nmask                                               # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_8-10_2025/Zm_nmask.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_8-10_2025/Zm_nmask.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_8-10_2025
ml BEDTools/2.31.1-GCC-13.3.0
bedtools maskfasta -fi Zm_B73_trim.fa -fo Zm_B73_A188_N-masked_genome.fa -bed syri.vcf -fullHeader