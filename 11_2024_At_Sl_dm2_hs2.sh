#!/bin/bash
#SBATCH --job-name=At_Sl_demultiplex2_hisat                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=6:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_At/At_Sl_dm2_hs2.out                  # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_At/At_Sl_dm2_hs2.err                   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_At/
mkdir Demultiplexed2
ml Miniconda3/23.5.2-0
source activate /home/jms53460/Fastq-Multx

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    if [ ! -f "Demultiplexed2/""$file2""_1s.fastq.gz" ]; then
        module load fastp/0.23.2-GCC-11.3.0
	    fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed2/umi_""$file2""_R1.fastq.gz" -O "Demultiplexed2/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

	    fastq-multx -b -B "CELSeq_barcodes.txt" -m 0 "Demultiplexed2/umi_""$file2""_R2.fastq.gz" "Demultiplexed2/umi_""$file2""_R1.fastq.gz" -o "Demultiplexed2/""$file2""_%_R2.fastq.gz" "Demultiplexed2/""$file2""_%.fastq.gz"  # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location

	    find "Demultiplexed2/" -name "umi_*" -delete
	    find "Demultiplexed2/" -name "*_R2*" -delete
    fi
done
conda deactivate

ml HISAT2/3n-20201216-gompi-2022a
mkdir hisat2_out2

for file in "Demultiplexed2/"*.fastq*
do
	file2="${file:15:-9}"

if [ ! -f "hisat2_out2/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.3.0
	fastp -w 6 -i "$file" -o "hisat2_out2/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

ml SAMtools/1.16.1-GCC-11.3.0
for file in "hisat2_out2/"*s.fastq*
do
	file2="${file:12:-9}"

if [ ! -f "hisat2_out2/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x Ler_0_N-masked/merged_hisat2_index -U "hisat2_out2/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out2/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out2/""$file2""_unsorted.bam" -o "hisat2_out2/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out2/""$file2""_s.bam"
	
fi
done