#!/bin/bash
#SBATCH --job-name=Zm_SNPsplit                                            # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_SNPsplit.out         # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_SNPsplit.err          # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1

SNPsplit --conflicting -o Zm_SNPsplit --snp_file Zm_SNPs.tab hisat2_out/merged_s.bam