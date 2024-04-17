#!/bin/bash
#SBATCH --job-name=Ns_SNPsplit_genome_prep                         # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=50gb                                                 # Total memory for job
#SBATCH --time=6:00:00                                             # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_SNPsplit_genome_prep2.out        # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_SNPsplit_genome_prep2.err         # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1

SNPsplit_genome_preparation --vcf_file Ns.vcf.gz --reference_genome Ns_genome.fna_dir --strain merged_s.bam --skip_filtering
conda deactivate
