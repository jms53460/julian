#!/bin/bash
#SBATCH --job-name=At_demultiplex                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7-8_2024_At/At_dm3.out                  # Location of standard output file
#SBATCH --error=/scratch/jms53460/7-8_2024_At/At_dm3.err                   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7-8_2024_At
module load fastp/0.23.2-GCC-11.3.0

for file in Demultiplexed2/*_R1_*.gz; do
    file2="${file:0:-15}"

        module load fastp/0.23.2-GCC-11.3.0
        fastp -w 6 -i $file -I "$file2""R2_001.fastq.gz" -o "$file2""_R1.fastq.gz" -O "$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

done