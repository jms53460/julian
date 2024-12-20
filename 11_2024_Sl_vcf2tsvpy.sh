#!/bin/bash
#SBATCH --job-name=Sl_minimap_syri                                         # Job name
#SBATCH --partition=batch                                                  # Partition (queue) name
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=6                                                  # Number of cores per task
#SBATCH --mem=350gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_Sl/Sl_minimap_syri.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_Sl/Sl_minimap_syri.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_Sl

ml minimap2/2.28-GCCcore-12.3.0
minimap2 -t 6 -ax asm10 --eqx Sl_genome_12.fa Sp_genome_12.fa > Sp_aligned_to_Sl_trim.sam

ml SyRI/1.6.3
syri -c Sp_aligned_to_Sl_trim.sam -r Sl_genome_12.fa -q Sp_genome_12.fa -k -F S
