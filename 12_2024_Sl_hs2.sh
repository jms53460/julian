#!/bin/bash
#SBATCH --job-name=Hisat2                                                 # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=70gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/12_2024_Sl/Hs2.out                     # Location of standard output file
#SBATCH --error=/scratch/jms53460/12_2024_Sl/Hs2.err                      # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/12_2024_Sl

#module load fastp/0.23.2-GCC-11.3.0
mkdir hisat2_out2
cp hisat2_out/*fastq.gz hisat2_out2
cp /scratch/jms53460/12_18_2024_Sl/*fastq.gz hisat2_out2
#for file in Demultiplexed/*s.fastq.gz; do
#	file2="${file:14:-9}"

#if [ ! -f "hisat2_out/""$file2"".bam" ]; then

#	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

#fi
#done

ml SAMtools/1.16.1-GCC-11.3.0
for file in hisat2_out2/*s.fastq.gz
do
	file2="${file:12:-9}"

if [ ! -f "hisat2_out2/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x Sl_N-masked_genome2_index -U "hisat2_out2/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out2/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out2/""$file2""_unsorted.bam" -o "hisat2_out2/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out2/""$file2""_s.bam"
	
fi
done