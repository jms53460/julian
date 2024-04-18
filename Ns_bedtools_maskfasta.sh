#!/bin/bash
#SBATCH --job-name=Ns_bedtools_maskfasta                           # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=1                                          # Number of cores per task
#SBATCH --mem=50gb                                                 # Total memory for job
#SBATCH --time=6:00:00                                             # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_bedtools_maskfasta.out          # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_bedtools_maskfasta.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460
ml BEDTools/2.30.0-GCC-12.2.0
bedtools maskfasta -fi Ns_genome.fna -fo Ns_N-masked_genome.fa -bed Ns.vcf.gz -fullHeader
