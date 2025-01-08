#!/bin/bash
#SBATCH --job-name=Ns_DNA_vcf                                      # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=50gb                                                 # Total memory for job
#SBATCH --time=12:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_vcf.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_vcf.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
module load SAMtools/1.16.1-GCC-11.3.0
samtools index -@ 6 hisat2_out/TW136_s.bam
samtools index -@ 6 hisat2_out/TW137_s.bam

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f Ns_genome.fna hisat2_out/TW136_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW136.vcf.gz
bcftools index TW136.vcf.gz

bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f Ns_genome.fna hisat2_out/TW137_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW137.vcf.gz
bcftools index TW137.vcf.gz
