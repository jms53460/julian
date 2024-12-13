#!/bin/bash
#SBATCH --job-name=At_Sl_demultiplex                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=48:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_At/At_Sl_dm3.out                  # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_At/At_Sl_dm3.err                   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_At/
mkdir Demultiplexed3
ml Miniconda3/23.5.2-0
source activate /home/jms53460/demultiplex

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    if [ ! -f "Demultiplexed3/""$file2""_1s.fastq.gz" ]; then
      
	    demultiplex demux -r -e 16 -m 0 -p Demultiplexed3 "CELSeq_barcodes.txt" "$file" "Raw_Data/""$file2""_R2_001.fastq.gz" # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location

    fi
done
conda deactivate