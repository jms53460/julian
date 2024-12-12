#!/bin/bash
#SBATCH --job-name=At_Sl_demultiplex                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_At/At_Sl_dm.out                  # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_At/At_Sl_dm.err                   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_At/
mkdir Demultiplexed
ml Miniconda3/23.5.2-0
source activate /home/jms53460/Fastq-Multx

module load fastp/0.23.2-GCC-11.3.0
if [ ! -f "Demultiplexed/""$file2""_1s.fastq.gz" ]; then
    for file in Raw_Data/*_R1_*.gz; do
        filename=$(basename "$file")
        file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

	    fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed/fastp_""$file2""_R1.fastq.gz" -O "Demultiplexed/fastp_""$file2""_R2.fastq.gz" -A -Q -L -G

	    fastq-multx -b -B "CELSeq_barcodes.txt" -m 0 "Demultiplexed/fastp_""$file2""_R2.fastq.gz" "Demultiplexed/fastp_""$file2""_R1.fastq.gz" -o "Demultiplexed/""$file2""_%_R2_dm.fastq.gz" "Demultiplexed/""$file2""_%_R1_dm.fastq.gz"  # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location
    done

    for file in Demultiplexed/*dm.fastq.gz; do
        file2="${file:0:-15}
        fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_dm.fastq.gz" -o "Demultiplexed/""$file2"".fastq.gz" -O "Demultiplexed/""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI
        
        find "Demultiplexed/" -name "*_dm*" -delete
        find "Demultiplexed/" -name "fastp_*" -delete
    done
fi

conda deactivate
