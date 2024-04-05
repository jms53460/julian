#!/bin/bash
#SBATCH --job-name=Ns_fastp_hisat2                      # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_fh2.out              # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_fh2.err               # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460
ml HISAT2/3n-20201216-gompi-2022a
hisat2-build Ns_genome.fna Ns_hisat2_index
mkdir hisat2_out

cd /home/jms53460/Mapped_data
for file in "demultiplexed/"*.fastq*
do
	file2="${file:14:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.3.0
	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 6 --dta -x /home/jms53460/Ns_hisat2_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 8 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2"".bam"
	
fi
done