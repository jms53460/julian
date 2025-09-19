#!/bin/bash
#SBATCH --job-name=Dm                                                     # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_8_2025/Dm.out                     # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_8_2025/Dm.err                      # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_8_2025/
mkdir Demultiplexed
ml Mamba/23.11.0-0
source activate /home/jms53460/Fastq-Multx

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    if [ ! -f "Demultiplexed/""$file2""_1s.fastq.gz" ]; then
        module load fastp/0.23.4-GCC-13.2.0
	    #Move UMI to header
        fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed/umi_""$file2""_R1.fastq.gz" -O "Demultiplexed/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI
        
        #Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location
	    fastq-multx -b -B "CELSeq_barcodes.txt" -m 0 "Demultiplexed/umi_""$file2""_R2.fastq.gz" "Demultiplexed/umi_""$file2""_R1.fastq.gz" "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed/""$file2""_%_R2.fastq.gz" "Demultiplexed/""$file2""_%.fastq.gz" "Demultiplexed/""$file2""_%_umi.fastq.gz"

    fi
done
conda deactivate

mkdir SRA_upload
module load fastp/0.23.4-GCC-13.2.0

for file in Demultiplexed/*s.fastq.gz; do
	file2="${file:14:-9}"

    #Trim UMI containing read to only contain the UMI
    fastp -w 6 -B 10 -i "Demultiplexed/""$file2"".fastq.gz" -I "Demultiplexed/""$file2""_umi.fastq.gz" -o "SRA_upload/""$file2"".fastq.gz" -O "SRA_upload/""$file2""_umi.fastq.gz" -A -Q -L -G
done