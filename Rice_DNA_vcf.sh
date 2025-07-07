#!/bin/bash
#SBATCH --job-name=Rice_DNA_vcf                                      # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=200gb                                                # Total memory for job
#SBATCH --time=24:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_DNA/Rice_vcf.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_DNA/Rice_vcf.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_DNA
module load SAMtools/1.16.1-GCC-11.3.0
samtools index -@ 6 hisat2_out/NRE1-2_s.bam

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f Nipponbare.fna hisat2_out/NRE1-2_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > NRE1-2.vcf.gz
bcftools index NRE1-2.vcf.gz