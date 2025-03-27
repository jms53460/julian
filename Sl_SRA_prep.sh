#!/bin/bash
#SBATCH --job-name=SRA_prep                                               # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/3_2025_Sl/SRA_prep.out                 # Location of standard output file
#SBATCH --error=/scratch/jms53460/3_2025_Sl/SRA_prep.err                  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/3_2025_Sl/
mkdir SRA_upload
module load fastp/0.23.2-GCC-11.3.0

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    #Trim UMI containing read to only contain the UMI
    fastp -w 6 -B 10 -i "Demultiplexed/""$file2""_%.fastq.gz" -I "Demultiplexed/""$file2""_%_umi.fastq.gz" -o "SRA_upload/""$file2""_%.fastq.gz" -O "SRA_upload/""$file2""_%_umi.fastq.gz"
done