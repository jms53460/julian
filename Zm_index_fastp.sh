#!/bin/bash
#SBATCH --job-name=Zm_index_fastp                                         # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_index_fastp.out      # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_index_fastp.err       # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022
ml HISAT2/3n-20201216-gompi-2022a
hisat2-build B73_v5_genome.fa B73_v5_hisat2_index
mkdir hisat2_out

for file in "Mapped_Data/demultiplexed2/"*.fastq*
do
	file2="${file:27:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.3.0
	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA
    
fi
done