#!/bin/bash
#SBATCH --job-name=Features_UMIs                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=1                                                     # Number of cores per task
#SBATCH --mem=250gb                                                            # Total memory for job
#SBATCH --time=6:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/12_2024_Sl/Features_UMIs.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/12_2024_Sl/Features_UMIs.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/12_2024_Sl
mkdir UMIcounts4
module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4

umi_tools count --per-contig -I "bams2/S14-25_S189_L007_10s_SNPsplit" -S "UMIcounts4/S14-25_S189_L007_10s_SNPsplit.tsv"
