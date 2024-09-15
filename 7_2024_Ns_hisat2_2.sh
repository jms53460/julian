#!/bin/bash
#SBATCH --job-name=Ns_hisat2_2                                              # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=2:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7_2024_Ns/Ns_hs2_2.out                 # Location of standard output file
#SBATCH --error=/scratch/jms53460/7_2024_Ns/Ns_hs2_2.err                  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7_2024_Ns
mkdir hisat2_out

for file in "Demultiplexed/"*.fastq*
do
	file2="${file:14:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.3.0
	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

mkdir hisat2_out2
cp hisat2_out/*s.fastq.gz hisat2_out2
ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0
for file in "hisat2_out2/"*s.fastq*
do
	file2="${file:12:-9}"

if [ ! -f "hisat2_out2/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x Ns_hisat2_index -U "hisat2_out2/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out2/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out2/""$file2""_unsorted.bam" -o "hisat2_out2/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out2/""$file2""_s.bam"
	
fi
done
