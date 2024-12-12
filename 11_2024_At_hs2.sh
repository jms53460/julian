#!/bin/bash
#SBATCH --job-name=At_hisat2                                              # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=2:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_At/At_hs2.out                 # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_At/At_hs2.err                  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_At/
mkdir hisat2_out
ml fastp/0.23.2-GCC-11.3.0
ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0

for file in "Demultiplexed/"*s__R1.fastq.gz
do
	file2="${file:14:-13}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA
    hisat2 -p 6 --dta -x Ler_0_N-masked/merged_hisat2_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out/""$file2""_s.bam"
	
fi
done