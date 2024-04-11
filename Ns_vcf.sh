#!/bin/bash
#SBATCH --job-name=Ns_vcf_snp                           # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_vcf_snp.out          # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_vcf_snp.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

module load SAMtools/1.16.1-GCC-11.3.0
samtools index -@ 6 hisat2_out/merged_s.bam

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 --min-MQ 60 -f Ns_genome.fna merged_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > Ns.vcf.gz
bcftools index Ns.vcf.gz