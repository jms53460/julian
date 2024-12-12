#!/bin/bash
#SBATCH --job-name=At_Sl_demultiplex                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                       # Total memory for job
#SBATCH --time=2:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_At/At_Sl_dm3.out                  # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_At/At_Sl_dm3.err                   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_At/
mkdir Demultiplexed3
ml Miniconda3/23.5.2-0
source activate /home/jms53460/Fastq-Multx

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    if [ ! -f "Demultiplexed3/""$file2""_1s.fastq.gz" ]; then
        module load fastp/0.23.2-GCC-11.3.0
	    fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed3/umi_""$file2""_R1.fastq.gz" -O "Demultiplexed3/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

	    fastq-multx -b -B "CELSeq_barcodes.txt" -m 0 "Demultiplexed3/umi_""$file2""_R2.fastq.gz" "Demultiplexed3/umi_""$file2""_R1.fastq.gz" -o "Demultiplexed3/""$file2""_%_R2.fastq.gz" "Demultiplexed3/""$file2""_%.fastq.gz"  # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location

	    find "Demultiplexed3/" -name "umi_*" -delete

    fi
done
conda deactivate
