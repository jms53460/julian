#!/bin/bash
#SBATCH --job-name=Sl_samtools_bcftools                                    # Job name
#SBATCH --partition=batch                                                  # Partition (queue) name
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=6                                                  # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_Sl/Sl_samtools_bcftools.out     # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_Sl/Sl_samtools_bcftools.err      # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_Sl

ml SAMtools/1.18-GCC-12.3.0
samtools view -bS Sp_aligned_to_Sl_trim.sam > Sp_aligned_to_Sl_trim.bam
samtools sort Sp_aligned_to_Sl.bam 

ml BCFtools/1.18-GCC-12.3.0
bcftools mpileup -Ou --threads 6 --min-MQ 60 -f Sl_genome_12.fa Sp_aligned_to_Sl.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > Sl.vcf.gz
bcftools index Sl.vcf.gz
