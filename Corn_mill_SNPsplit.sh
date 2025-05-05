#!/bin/bash
#SBATCH --job-name=Corn_mill_SNPsplit
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=8                                                  # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Corn_mill/Corn_mill_SNP.out             # Location of standard output file
#SBATCH --error=/scratch/jms53460/Corn_mill/Corn_mill_SNP.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Corn_mill

mkdir SNPsplit
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "mapped/"*.bam
do
    file2="${file:7:-4}"

    SNPsplit --conflicting -o SNPsplit --snp_file Zm_W22_w_B73_SNPs.tab "$file"
    samtools sort -@ 6 SNPsplit/"$file2".allele_flagged.bam -o SNPsplit/"$file2"SNPsplit.bam
    
done