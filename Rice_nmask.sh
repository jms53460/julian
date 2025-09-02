#!/bin/bash
#SBATCH --job-name=Rice_nmask                                             # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_DNA2/Rice_nmask.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_DNA2/Rice_nmask.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_DNA2
ml BEDTools/2.31.1-GCC-13.3.0
for file in "Raw_Data/"*_R1_001.fastq.gz
do
	file2="${file:9:-16}"
bedtools maskfasta -fi Nipponbare_12.fna -fo "$file2"_N-masked_12.fa -bed "$file2"_snps.vcf -fullHeader
done