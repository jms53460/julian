#!/bin/bash
#SBATCH --job-name=RiceRemap_nmask                                             # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/RiceRemap/Rice_nmask.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/RiceRemap/Rice_nmask.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/RiceRemap
ml BEDTools/2.31.1-GCC-13.3.0
for file in "Raw_Data_DNA/"*_R1_001.fastq.gz
do
	file2="${file:13:-16}"
bedtools maskfasta -fi Rice_r7.fa -fo "$file2"_N-masked_Remap.fa -bed "$file2"_snps.vcf -fullHeader
done