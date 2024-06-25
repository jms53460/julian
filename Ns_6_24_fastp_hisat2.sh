#!/bin/bash
#SBATCH --job-name=Ns_6_24_fastp_hisat2                               # Job name
#SBATCH --partition=batch                                             # Partition (queue) name
#SBATCH --ntasks=1                                                    # Single task job
#SBATCH --cpus-per-task=6                                             # Number of cores per task
#SBATCH --mem=50gb                                                    # Total memory for job
#SBATCH --time=12:00:00                                               # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/June2024SeqNs_fh2.out              # Location of standard output file
#SBATCH --error=/scratch/jms53460/June2024Seq/Ns_fh2.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                  # Where to send mail
#SBATCH --mail-type=END,FAIL                                          # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/June2024Seq

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

        module load fastp/0.23.2-GCC-11.3.0
	    fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Mapped_Data/umi_""$file2""_R1.fastq.gz" -O "Mapped_Data/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

	    find "Mapped_Data/" -name "*_R2*" -delete
done

ml HISAT2/3n-20201216-gompi-2022a

for file in "Mapped_Data/"*.fastq*
do
	file2="${file:12:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.3.0
	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 6 --dta -x Ns_hisat2_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
	
fi
done