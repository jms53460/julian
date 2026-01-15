#!/bin/bash
#SBATCH --job-name=RiceRemap_fastp_hisat2                             # Job name
#SBATCH --partition=batch                                             # Partition (queue) name
#SBATCH --ntasks=1                                                    # Single task job
#SBATCH --cpus-per-task=6                                             # Number of cores per task
#SBATCH --mem=200gb                                                   # Total memory for job
#SBATCH --time=12:00:00                                               # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/RiceRemap/Rice_fh2.out             # Location of standard output file
#SBATCH --error=/scratch/jms53460/RiceRemap/Rice_fh2.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                  # Where to send mail
#SBATCH --mail-type=END,FAIL                                          # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/RiceRemap

ml HISAT2/2.2.1-gompi-2023a
hisat2-build Rice_r7.fa Rice_r7_index

mkdir hisat2_out_DNA
module load fastp/0.23.4-GCC-13.2.0
module load SAMtools/1.21-GCC-13.3.0
for file in "Raw_Data_DNA/"*_R1_001.fastq.gz
do
	file2="${file:13:-16}"

	fastp -w 6 -i "$file" -I "Raw_Data_DNA/""$file2""_R2_001.fastq.gz" -o "hisat2_out_DNA/""$file2""_R1_001.fastq.gz" -O "hisat2_out_DNA/""$file2""_R2_001.fastq.gz"

	hisat2 -p 6 --dta -x Rice_r7_index -1 "hisat2_out_DNA/""$file2""_R1_001.fastq.gz" -2 "hisat2_out_DNA/""$file2""_R2_001.fastq.gz" | samtools view -bS -> "hisat2_out_DNA/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out_DNA/""$file2""_unsorted.bam" -o "hisat2_out_DNA/""$file2""_s.bam"
done