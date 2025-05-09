#!/bin/bash
#SBATCH --job-name=Corn_mill_SNPsplit
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=8                                                  # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Corn_mill/Corn_mill_SNP2.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/Corn_mill/Corn_mill_SNP2.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Corn_mill

mkdir SNPsplit2
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "mapped/"*_.bam
do
    file2="${file:7:-5}"
    samtools sort -@ 8 mapped/"$file2""_.bam" -o mapped/"$file2""_s.bam"
    samtools index -@ 8 mapped/"$file2""_s.bam"
    SNPsplit --conflicting -o SNPsplit2 --snp_file Zm_W22_w_B73_SNPs.tab mapped/"$file2""_s.bam"
    samtools sort -@ 6 SNPsplit2/"$file2"_s.allele_flagged.bam -o SNPsplit2/"$file2"SNPsplit.bam
    
done