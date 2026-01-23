#!/bin/bash
#SBATCH --job-name=TomatoRemap_minimap_syri_vcf2tsv                        # Job name
#SBATCH --partition=batch                                                  # Partition (queue) name
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=6                                                  # Number of cores per task
#SBATCH --mem=350gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/TomatoRemap/Sl_minimap_syri.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/TomatoRemap/Sl_minimap_syri.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/TomatoRemap

ml minimap2/2.28-GCCcore-13.2.0
minimap2 -t 6 -ax asm10 --eqx Sl_v4_12.fa Sp_genome_12.fa > Sp_aligned_to_Sl_trim.sam

ml SyRI/1.7.1-foss-2023a
syri -c Sp_aligned_to_Sl_trim.sam -r Sl_v4_12.fa -q Sp_genome_12.fa -k -F S

ml Mamba/23.11.0-0
source activate /home/jms53460/vcf2tsvpy
vcf2tsvpy --input_vcf syri.vcf --out_tsv Sl_vcf_table.tsv 
conda deactivate