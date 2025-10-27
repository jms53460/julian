#!/bin/bash
#SBATCH --job-name=Ran_features_rev                                       # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=300gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Ran/Ran_features_rev.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/Ran/Ran_features_rev.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Ran
mkdir featurecounts_rev

ml Mamba/23.11.0-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 2 -a JEC21_56.gff -t 'protein_coding_gene' -g 'ID' -o featurecounts_rev/read_counts.tab --readExtension5 500 -R BAM Mapped/*2*.bam

conda deactivate
