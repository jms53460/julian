#!/bin/bash
#SBATCH --job-name=Rice_DNA_rmdup                                  # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=100gb                                                # Total memory for job
#SBATCH --time=12:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_DNA/Rice_rmdup.out         # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_DNA/Rice_rmdup.err          # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_DNA
module load SAMtools/1.16.1-GCC-11.3.0
mkdir rmdup

samtools rmdup hisat2_out/NRE1-2_s.bam rmdup/NRE1-2_rmdup.bam
samtools index -@ 6 rmdup/NRE1-2_rmdup.bam