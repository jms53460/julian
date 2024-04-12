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
bcftools mpileup -Ou --threads 6 -g 10 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD -E -Q 0 -pm 3 -F 0.25 -d 500 -f Ns_genome.fna merged_s.bam | bcftools call -Ou -a GQ,GP -mAv -p 0.99 --threads 6 | bcftools norm -Ou --threads 6 -f Ns_genome.fna -m +indels | bcftools filter -Oz -e 'QUAL<40 || DP<10' --threads 6 > Ns.3_vcf.gz
bcftools view -Ou --threads 6 -i 'FORMAT/FI[*] = 1' Ns.3_vcf.gz | bcftools annotate -Ou --threads 6 -x INFO/VDB,INFO/SGB,INFO/RPBZ,INFO/MQBZ,INFO/MQBZ,INFO/MQSBZ,INFO/BQBZ,INFO/SCBZ,INFO/FS,INFO/MQOF,INFO/AC,INFO/AN,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,FORMAT/GP | bcftools view --threads 6 -a -Oz -o Ns.4_vcf.gz
bcftools annotate -x INFO/MQ0F -Oz --threads 6-o Ns.4_vcf.gz Ns.3_vcf.gz
bcftools index Ns.4_vcf.gz