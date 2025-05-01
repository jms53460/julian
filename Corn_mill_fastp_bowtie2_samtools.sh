#!/bin/bash
#SBATCH --job-name=JuneSequencing2024W22
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=8                                                  # Number of cores per task
#SBATCH --mem=200gb                                                        # Total memory for job
#SBATCH --time=65:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Corn_mill/Corn_mill_map.out             # Location of standard output file
#SBATCH --error=/scratch/jms53460/Corn_mill/Corn_mill_map.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Corn_mill
ml fastp/0.23.4-GCC-12.3.0
ml Bowtie2/2.5.2-GCC-11.3.0
ml SAMtools/1.18-GCC-12.3.0
ml UMI-tools/1.1.2-foss-2022a-Python-3.10.4

bowtie2-build Zm_W22_w_B73_N-masked.fa Zm_W22_w_B73_N-masked

for file in rawdata/*R1*
do
        file2="${file:8:-15}"
        fastp -w 8 -i $file -o processed/"$file2""R1A.fastq.gz" -f 23 -A -G -Q -L
        fastp -i processed/"$file2""R1A.fastq.gz" -I rawdata/"$file2""R2_001.fastq.gz" -o processed/"$file2""R1B.fastq.gz" -O processed/"$file2""R2B.fastq.gz" -U --umi_loc read1 --umi_len 6 --umi_prefix TIR -A -G -Q -L
        fastp -i processed/"$file2""R1B.fastq.gz" -I processed/"$file2""R2B.fastq.gz" -o processed/"$file2""R1C.fastq.gz" -O processed/"$file2""R2C.fastq.gz" -U --umi_loc read2 --umi_len 8 --umi_prefix UMI -A -G -Q -L
        fastp -i processed/"$file2""R1C.fastq.gz" -I processed/"$file2""R2C.fastq.gz" -o filtered/"$file2""R1.fastq.gz" -O filtered/"$file2""R2.fastq.gz" -U --umi_loc read2 --umi_len 11 --umi_prefix BC --length_required 40 --trim_poly_x --cut_tail --adapter_fasta /scratch/jms53460/Corn_mill/adapter2.1.fa
bowtie2 --threads 8 -x Zm_W22_w_B73_N-masked --phred33 -X 1000 --no-mixed --no-discordant -1 filtered/"$file2""R1.fastq.gz" -2 filtered/"$file2""R2.fastq.gz" | samtools view -@ 8 -b -o mapped/"$file2"".bam"
        samtools sort -@ 8 mapped/"$file2"".bam" -o mapped/"$file2"".bam"
        samtools index -@ 8 mapped/"$file2"".bam"
done
