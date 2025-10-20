#!/bin/bash
#SBATCH --job-name=Zm_minimap_syri                                         # Job name
#SBATCH --partition=batch                                                  # Partition (queue) name
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=6                                                  # Number of cores per task
#SBATCH --mem=350gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_8-10_2025/Zm_minimap_syri.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_8-10_2025/Zm_minimap_syri.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_8-10_2025

ml minimap2/2.28-GCCcore-13.2.0
minimap2 -t 6 -ax asm20 --eqx Zm-B73-REFERENCE-NAM-5.0.fa Zm_A188.fa > A188_aligned_to_B73.sam

ml SyRI/1.7.1-foss-2023a
syri -c A188_aligned_to_B73.sam -r Zm-B73-REFERENCE-NAM-5.0.fa -q Zm_A188.fa -k -F S