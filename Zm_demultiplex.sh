#!/bin/bash
#SBATCH --job-name=Zm_demultiplex                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=200gb                                                       # Total memory for job
#SBATCH --time=48:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_dm.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_dm.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022
ml Miniconda3/23.5.2-0
source activate /scratch/jms53460/Fastq-Multx/

for file in Raw_Data/*_R1*.gz; do
    file2="${file:9:-12}"

    if [ ! -f "Mapped_Data/demultiplexed/""$file2""_dT-1s.fastq.gz" ]; then
        module load fastp/0.23.2-GCC-11.3.0
	    fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2.fastq.gz" -o "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" -O "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" -A -Q -L --umi --umi_loc read1 --umi_len 10 --umi_prefix UMI

	    fastq-multx -b -B "CELSeq_barcodes.txt" -m 0 "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" -o "Mapped_Data/demultiplexed/""$file2""_%_R1.fastq.gz" "Mapped_Data/demultiplexed/""$file2""_%.fastq.gz"  # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location

	    find "Mapped_Data/demultiplexed/" -name "umi_*" -delete
	    find "Mapped_Data/demultiplexed/" -name "*_R1*" -delete
    fi
done
conda deactivate