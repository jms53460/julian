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
mkdir Demultiplexed
ml Miniconda3/23.5.2-0
source activate /home/jms53460/Fastq-Multx

for file in Raw_Data/*R1*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

	    fastq-multx -B "CELSeq_barcodes.txt" -m 0 "Raw_Data/""$file2""_R2_001.fastq.gz" "Raw_Data/""$file2""_R1_001.fastq.gz" -o "Demultiplexed/""$file2""_%_R2_dm.fastq.gz" "Demultiplexed/""$file2""_%_R1_dm.fastq.gz"  # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location
	
done


module load fastp/0.23.2-GCC-11.3.0
for file in Demultiplexed/*_R1_*.gz; do
    file2="${file:0:-14}"

        fastp -w 6 -i $file -I "$file2""R2_dm.fastq.gz" -o "$file2""R1.fastq.gz" -O "$file2""R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

done
conda deactivate
