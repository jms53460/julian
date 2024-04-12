#!/bin/bash
#SBATCH --job-name=Ns_vcf3                              # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_vcf3.out             # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_vcf3.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460
module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 --min-MQ 60 -f Ns_genome.fna merged_s.bam | bcftools call -Ou -m -v -f GQ,GP --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' --threads 6 > Ns.3_vcf.gz
bcftools index Ns.3_vcf.gz