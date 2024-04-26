#!/bin/bash
#SBATCH --job-name=Zm_download_v4                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=3:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_download4.out         # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_download4.err          # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022

curl -s https://stacks.stanford.edu/file/druid:js115sm3463/4o-A188_B73.AB.final.vcf > v4_A188_B73.vcf
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/A188_SNPs.tab.gz > v4_SNPs.tab.gz
curl -s https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz > v4_B73.fa
