#!/bin/bash
#SBATCH --job-name=At_hisat2                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_At/At_hs2.out                  # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_At/At_hs2.err                   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_At/

module load fastp/0.23.2-GCC-11.3.0
mkdir UMI_moved
for file in Demultiplexed3/A254-266*R1*s.fastq.gz; do
        file2="${file:39:-9}"
        fastp -w 6 -i "$file" -I "Demultiplexed3/A254-266_S2_L002_R2_001_""$file2"".fastq.gz" -o "UMI_moved/A254-266_"$file2""_R1.fastq.gz" -O "UMI_moved/A254-266_"$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI
    done
for file in Demultiplexed3/S1-8_A267-277*R1*s.fastq.gz; do
        file2="${file:44:-9}"
        fastp -w 6 -i "$file" -I "Demultiplexed3/S1-8_A267-277_S3_L002_R2_001_""$file2"".fastq.gz" -o "UMI_moved/S1-8_A267-277_"$file2""_R1.fastq.gz" -O "UMI_moved/S1-8_A267-277_"$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI
    done

mkdir hisat2_out

for file in "UMI_moved/"*_R1.fastq.gz
do
	file2="${file:10:-12}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0
for file in "hisat2_out/"*s.fastq.gz
do
	file2="${file:11:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x Ler_0_N-masked/merged_hisat2_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out/""$file2""_s.bam"
	
fi
done