#!/bin/bash
#SBATCH --job-name=Zm_mpileup                                             # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=12                                                # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_mpileup.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_mpileup.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022
module load SAMtools/1.16.1-GCC-11.3.0
samtools view -@ 12 -bS A188_aligned_to_B73.sam > A188_aligned_to_B73.bam
samtools index -@ 12 A188_aligned_to_B73.bam

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 12 --min-MQ 60 -f B73_v5_genome.fa A188_aligned_to_B73.bam | bcftools call -Ou -m -v --threads 12 | bcftools filter --threads 12 -Oz -e 'QUAL<40 || DP<10' > Zm_mpileup.vcf.gz
bcftools index --threads 12 Zm_mpileup.vcf.gz