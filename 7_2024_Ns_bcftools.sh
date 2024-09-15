#!/bin/bash
#SBATCH --job-name=Ns_bcftools                                              # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=2:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7_2024_Ns/Ns_bcftools.out                 # Location of standard output file
#SBATCH --error=/scratch/jms53460/7_2024_Ns/Ns_bcftools.err                  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7_2024_Ns
module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 -d 1000 --min-MQ 60 --skip-indels -f Ns_genome.fna -b hisat2_out2/*_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > Ns_2.vcf.gz
bcftools index Ns_2.vcf.gz