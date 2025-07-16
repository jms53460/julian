#!/bin/bash
#SBATCH --job-name=SRA_prep                                               # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_7_2025/SRA_prep.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_7_2025/SRA_prep.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_7_2025/
mkdir SRA_upload
module load fastp/0.23.2-GCC-11.3.0

for file in Demultiplexed/*s.fastq.gz; do
	file2="${file:14:-9}"

    #Trim UMI containing read to only contain the UMI
    fastp -w 6 -B 10 -i "Demultiplexed/""$file2"".fastq.gz" -I "Demultiplexed/""$file2""_umi.fastq.gz" -o "SRA_upload/""$file2"".fastq.gz" -O "SRA_upload/""$file2""_umi.fastq.gz" -A -Q -L -G
done