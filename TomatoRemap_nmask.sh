#!/bin/bash
#SBATCH --job-name=TomatoRemap_nmask                                      # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=6:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/TomatoRemap/Sl_nmask.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/TomatoRemap/Sl_nmask.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/TomatoRemap
ml BEDTools/2.31.1-GCC-13.3.0
bedtools maskfasta -fi Sl_v4_12.fa -fo Sl_v4_N-masked_genome.fa -bed syri.vcf -fullHeader

ml HISAT2/2.2.1-gompi-2023a
hisat2-build Sl_v4_N-masked_genome.fa Sl_v4_N-masked_genome_index