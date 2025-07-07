#!/bin/bash
#SBATCH --job-name=Rice_fastp_hisat2                                  # Job name
#SBATCH --partition=batch                                             # Partition (queue) name
#SBATCH --ntasks=1                                                    # Single task job
#SBATCH --cpus-per-task=6                                             # Number of cores per task
#SBATCH --mem=50gb                                                    # Total memory for job
#SBATCH --time=12:00:00                                               # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_DNA/Rice_fh2.out              # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_DNA/Rice_fh2.err               # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                  # Where to send mail
#SBATCH --mail-type=END,FAIL                                          # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_DNA

ml HISAT2/3n-20201216-gompi-2022a
hisat2-build Nipponbare.fna Nipponbare_index

mkdir hisat2_out
for file in "Raw_Data/"*_R1_001.fastq.gz
do
	file2="${file:9:-16}"

	module load fastp/0.23.2-GCC-11.3.0
	fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "hisat2_out/""$file2""_R1_001.fastq.gz" -O "hisat2_out/""$file2""_R2_001.fastq.gz"

	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 6 --dta -x Nipponbare_index -1 "hisat2_out/""$file2""_R1_001.fastq.gz" -2 "hisat2_out/""$file2""_R2_001.fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
done