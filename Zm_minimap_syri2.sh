#!/bin/bash
#SBATCH --job-name=Zm_minimap_syri                                         # Job name
#SBATCH --partition=batch                                                  # Partition (queue) name
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=6                                                  # Number of cores per task
#SBATCH --mem=350gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_minimap_syri2.out     # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_minimap_syri2.err      # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022

ml minimap2/2.26-GCCcore-12.2.0
minimap2 -t 6 -ax asm5 --eqx B73_v5_genome_trim.fa A188_genome_trim.fa > A188_aligned_to_B73_trim.sam

ml SyRI/1.6.3
syri -c A188_aligned_to_B73_trim.sam -r B73_v5_genome_trim.fa -q A188_genome_trim.fa -k -F S
