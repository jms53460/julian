#!/bin/bash
#SBATCH --job-name=SNPsplit                                                      # Job name
#SBATCH --partition=batch                                                           # Partition (queue) name
#SBATCH --ntasks=1                                                                  # Single task job
#SBATCH --cpus-per-task=6                                                           # Number of cores per task
#SBATCH --mem=50gb                                                                  # Total memory for job
#SBATCH --time=6:00:00                                                              # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_Sl/SNPsplit.out                       # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_Sl/SNPsplit.err                        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                                # Where to send mail
#SBATCH --mail-type=END,FAIL                                                        # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_Sl
mkdir SNPsplit2
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "hisat2_out2/"*_s.bam
do
    file2="${file:12:-6}"

    SNPsplit --conflicting -o SNPsplit2 --snp_file Sl_SNPs.tab "$file"
    samtools sort -@ 6 SNPsplit2/"$file2"_s.allele_flagged.bam -o SNPsplit2/"$file2"_SNPsplit.bam
    
done

for file in "SNPsplit2/"*_s.genome1.bam
do
    file2="${file:10:-14}"
    samtools sort -@ 6 "$file" -o SNPsplit2/"$file2"_SNPsplit_g1.bam
done

for file in "SNPsplit2/"*_s.genome2.bam
do
    file2="${file:10:-14}"
    samtools sort -@ 6 "$file" -o SNPsplit2/"$file2"_SNPsplit_g2.bam
done
